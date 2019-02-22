import numpy as np
import datetime
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText


class Counter(object):
    __slots__ = (
        "steps_per_temperature",
        "flips_per_step",
        "attempts_per_flip",
        "T_min",
        "T_max",
        "T_step",
        "T_counter",
    )

    def __init__(self, params):
        self.steps_per_temperature = int(params["stepsPerTemperature"])
        self.flips_per_step = int(params["flipsPerStep"])
        self.attempts_per_flip = int(params["attemptsPerFlip"])
        self.T_min = float(params["TMin"])
        self.T_max = float(params["TMax"])
        self.T_step = 0.1
        self.T_counter = 0

    def temperature(self):
        return self.T_min + self.T_counter * self.T_step

    def steps_iter(self):
        return range(self.steps_per_temperature)

    def flips_iter(self):
        return range(self.flips_per_step)

    def flip_attempts_iter(self):
        return range(self.attempts_per_flip)

    def update_temperature(self):
        self.T_counter += 1

    def should_continue(self):
        return self.T_min + self.T_counter * self.T_step <= self.T_max


class Measurements(object):
    __slots__ = (
        "sum_energy_sq",
        "sum_energy",
        "sum_magnetization_sq",
        "sum_magnetization",
    )

    def __init__(self):
        self.sum_energy_sq = 0
        self.sum_energy = 0
        self.sum_magnetization_sq = 0
        self.sum_magnetization = 0

    def zero(self):
        self.__init__()


class State(object):
    __slots__ = (
        "lattice",
        "counter",
        "lattice_size",
        "j_const",
        "k_const",
        "measurements",
    )

    def __init__(self, params):
        size = int(params["latticeSize"])
        self.lattice = generate_lattice(size)
        self.counter = Counter(params)
        self.lattice_size = size
        self.j_const = 1
        self.k_const = 1
        self.measurements = Measurements()


def generate_lattice(size):
    return np.random.choice([-1, 1], (size, size))


def random_cell(size):
    return np.random.randint(0, size, (2, 1))


def calc_d_energy(state, i, j):
    j_const = state.j_const
    k_const = state.k_const
    size = state.lattice_size
    lattice = state.lattice

    # d_energy = e_1 - e_0 =
    #  = ((-J) * (-spin) * neigbours) - ((-J) * spin * neighbours) =
    #  = -J * neighbours * ((-spin) - spin) =
    #  = -2 * -J * neighbours * spin =
    #  = 2 * J * spin * neighbours

    return (
        2
        * j_const
        * lattice[i, j]
        * (
            lattice[(i + 1 + size) % size, j]
            + lattice[(i - 1 + size) % size, j]
            + lattice[i, (j + 1 + size) % size]
            + lattice[i, (j - 1 + size) % size]
        )
    )


def is_flip_probable(state, i, j):
    d_energy = calc_d_energy(state, i, j)
    T = state.counter.temperature()
    k_const = state.k_const

    if d_energy < 0:
        return True
    else:
        rand = np.random.random()
        probability = np.exp(-d_energy / (T * k_const))

        return rand < probability


def flip_spin(state):
    lattice = state.lattice
    size = state.lattice_size

    for attempt in state.counter.flip_attempts_iter():
        i, j = random_cell(size)

        if is_flip_probable(state, i, j):
            lattice[i, j] *= -1

            break


def calc_energy(state):
    lattice = state.lattice

    return -j_const * np.sum(
        lattice * (np.roll(lattice, 1, axis=0) + np.roll(lattice, 1, axis=1))
    )


def calc_fluctuations(state):
    measurements = state.measurements
    n_steps = state.counter.steps_per_temperature
    T = state.counter.temperature()

    avg_energy_sq = measurements.sum_energy_sq / n_steps
    avg_energy = measurements.sum_energy / n_steps

    return (avg_energy_sq - avg_energy ** 2) / T


def calc_magnetization(state):
    measurements = state.measurements
    n_steps = state.counter.steps_per_temperature

    return measurements.sum_magnetization / n_steps


def calc_mag_susceptibility(state):
    measurements = state.measurements
    n_steps = state.counter.steps_per_temperature
    T = state.counter.temperature()

    avg_magnetization_sq = measurements.sum_magnetization_sq / n_steps
    avg_magnetization = measurements.sum_magnetization / n_steps

    return (avg_magnetization_sq - avg_magnetization ** 2) / T


def calculate(state):
    measurements = state.measurements

    fluctuations = calc_fluctuations(state)
    magnetization = calc_magnetization(state)
    mag_susceptibility = calc_mag_susceptibility(state)

    record = {
        "T": state.counter.temperature(),
        "sumEnergySq": measurements.sum_energy_sq,
        "sumEnergy": measurements.sum_energy,
        "fluctuations": fluctuations,
        "magnetization": magnetization,
        "magSusceptibility": mag_susceptibility,
    }

    # zero the measurements for the next temperature
    measurements.zero()

    return record


def measure_energy(state):
    lattice = state.lattice
    j_const = state.j_const

    return -j_const * np.sum(
        lattice * (np.roll(lattice, 1, axis=0) + np.roll(lattice, 1, axis=1))
    )


def measure_magnetization(state):
    return abs(np.mean(state.lattice))


def measure(state):
    measurements = state.measurements

    energy = measure_energy(state)
    magnetization = measure_magnetization(state)

    measurements.sum_energy_sq += energy ** 2
    measurements.sum_energy += energy
    measurements.sum_magnetization_sq += magnetization ** 2
    measurements.sum_magnetization += magnetization


def setup(params):
    state = State(params)
    flips_to_skip = int(params["flipsToSkip"])

    for flip in range(flips_to_skip):
        flip_spin(state)

    return state


def update(state):
    lattice = state.lattice
    counter = state.counter

    for step in counter.steps_iter():
        for flip in counter.flips_iter():
            flip_spin(state)

        measure(state)

    record = calculate(state)
    counter.update_temperature()

    return record, counter.should_continue()


def simulate(setup_func, params, update_func, measure_func):
    state = setup_func(params)
    results = []

    while True:
        record, should_continue = update_func(state)
        results.append(record)

        if not should_continue:
            break

    return results


def format_d_t(d_t):
    hours, r = divmod(d_t.seconds, 3600)
    minutes, seconds = divmod(r, 60)

    return "{:02d}h {:02d}m {:02d}s".format(hours, minutes, seconds)


def compose_results_file_path(prefix):
    date = datetime.datetime.now().strftime("%d-%m-%Y-%H:%M")
    filename = "{prefix}/results-{date}.txt".format(prefix=prefix, date=date)

    return filename


def format_record(record):
    return "{:>5.2f}{:>15.0f}{:>15.2f}{:>20.2f}{:>15.0f}{:>15.0f}".format(
        record["T"],
        record["fluctuations"],
        record["magnetization"],
        record["magSusceptibility"],
        record["sumEnergySq"],
        record["sumEnergy"],
    )


def compose_results(records, d_t):
    headers = "{:>5}{:>15}{:>15}{:>20}{:>15}{:>15}".format(
        "T",
        "fluctuations",
        "magnetization",
        "mag. susceptibility",
        "sum energy sq",
        "sum energy",
    )
    contents = "\n".join([format_record(record) for record in records])
    d_t = format_d_t(d_t)

    return "{headers}\n{contents}\n\ntime: {d_t}\n".format(
        headers=headers, contents=contents, d_t=d_t
    )


def save_results_to_file(results, prefix):
    with open(compose_results_file_path(prefix), "w+") as f:
        f.write(results)


def send_mail_notification(results, sender, receiver, password):
    msg = MIMEMultipart()
    msg["From"] = sender
    msg["To"] = receiver
    msg["Subject"] = "Simulation results"
    body = "Results:\n{}\n".format(results)
    msg.attach(MIMEText(body, "plain"))

    server = smtplib.SMTP("smtp.gmail.com", 587)
    server.starttls()
    server.login(sender, password)
    server.sendmail(sender, receiver, msg.as_string())
    server.quit()
