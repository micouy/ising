import numpy as np
import time
from time import time
import datetime


SIZE = 50
J = 1

MIN_TEMPERATURE = 9.0
MAX_TEMPERATURE = 10.0

MAX_FLIP_ATTEMPTS = 30
FLIPS_PER_UPDATE = 250  # around SIZE^2
UPDATES_PER_TEMPERATURE = 1000
UPDATES_TO_SKIP = 300

PRINT_RESULTS = False


def generate_lattice():
    return np.random.choice([-1, 1], (SIZE, SIZE))


def show_lattice(lattice):
    print "\n".join(
            ["".join(["##" if cell == -1 else " ." for cell in row]) for row in lattice]
        )


def random_cell():
    return (np.random.randint(0, SIZE), np.random.randint(0, SIZE))


def is_change_probable(d_energy, temperature):
    if d_energy < 0:
        return True
    else:
        rand = np.random.random()
        probability = np.exp(-d_energy / temperature)

        return rand < probability


def calc_d_energy(j_const, lattice, i, j):
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
            lattice[(i + 1 + SIZE) % SIZE, j]
            + lattice[(i - 1 + SIZE) % SIZE, j]
            + lattice[i, (j + 1 + SIZE) % SIZE]
            + lattice[i, (j - 1 + SIZE) % SIZE]
        )
    )


def flip_spin(lattice, i, j):
    lattice[i][j] *= -1


def calc_energy(j_const, lattice):
    return -j_const * np.sum(
        lattice * (np.roll(lattice, 1, axis=0) + np.roll(lattice, 1, axis=1))
    )


def flip(j_const, max_flip_attempts, temperature, lattice):
    for a in range(max_flip_attempts):
        i, j = random_cell()
        d_energy = calc_d_energy(j_const, lattice, i, j)

        if is_change_probable(d_energy, temperature):
            flip_spin(lattice, i, j)

            break


def update(j_const, max_flip_attempts, flips_per_update, temperature, lattice):
    for f in range(flips_per_update):
        flip(j_const, max_flip_attempts, temperature, lattice)


def temperature_steps(T_min, T_max):
    step = 0.1
    d_T = T_max - T_min
    n_steps = int(d_T / step)
    counter = 0

    while counter <= n_steps:
        yield T_min + step * counter

        counter += 1


def calc_avg_fluctuation_at_T(sum_energy_sq, sum_energy, updates_passed, T):
    avg_energy_sq = sum_energy_sq / updates_passed
    avg_energy = sum_energy / updates_passed

    return (avg_energy_sq - avg_energy ** 2) / T


def calc_magnetization(lattice):
    return np.mean(lattice)


def calc_avg_magnetization(sum_magnetization, updates_passed):
    return sum_magnetization / updates_passed


def display_latest_results(
    sum_energy_sq, sum_energy, sum_magnetization, updates_passed, T, t_start
):
    avg_fluctuation = calc_avg_fluctuation_at_T(
        sum_energy_sq, sum_energy, updates_passed, T
    )
    avg_magnetization = calc_avg_magnetization(sum_magnetization, updates_passed)

    t = time() - t_start
    print "Update #{}, T = {:.2f}, t = {:.2f}".format(updates_passed, T, t)
    print "Fluctuations: {:.0f}".format(avg_fluctuation)
    print "Magnetization: {:.2f}".format(avg_magnetization)


def simulate(temperature_params, cycle_params, j_const):
    T_min, T_max = temperature_params
    updates_to_skip, max_flip_attempts, flips_per_update, updates_per_temperature = (
        cycle_params
    )

    t = time()
    lattice = generate_lattice()
    results = {}

    for T in temperature_steps(T_min, T_max):
        # optinally we can generate a new lattice at each temperature step
        # lattice = generate_lattice()
        sum_energy_sq = 0
        sum_energy = 0
        sum_magnetization = 0

        for u in range(updates_to_skip):
            update(j_const, max_flip_attempts, flips_per_update, T, lattice)

        for u in range(updates_per_temperature):
            updates_passed = u + 1

            if updates_passed % 100 == 0 and PRINT_RESULTS:
                # show_lattice(lattice)
                display_latest_results(
                    sum_energy_sq, sum_energy, sum_magnetization, updates_passed, T, t
                )

            update(j_const, max_flip_attempts, flips_per_update, T, lattice)

            energy = calc_energy(j_const, lattice)
            magnetization = calc_magnetization(lattice)

            sum_energy_sq += energy ** 2
            sum_energy += energy
            sum_magnetization += magnetization

        results[T] = {
            "fluctuations": calc_avg_fluctuation_at_T(
                sum_energy_sq, sum_energy, updates_per_temperature, T
            ),
            "magnetization": calc_avg_magnetization(
                sum_magnetization, updates_per_temperature
            ),
        }

    return results


def compose_filename():
    date = datetime.datetime.now().strftime("%d-%m-%Y-%H:%M")
    filename = "results/results-{date}.txt".format(date=date)

    return filename


def format_record(record):
    T, data = record

    return "{:>15.2f}{:>20.0f}{:>15.2f}".format(T, data["fluctuations"], data["magnetization"])

def format_d_t(d_t):
    return "{:02d}:{:02d}".format(int(math.floor(d_t.seconds / 3600)), d_t.seconds % 3600)

def save_results(results, d_t):
    print results

    filename = compose_filename()

    with open(filename, "w+") as f:
        headers = "{:>15}{:>20}{:>15}".format("temperature", "fluctuations", "magnetization")
        records = map(format_record, results.items())
        contents = "\n".join(records)
        d_t = format_d_t(d_t)

        f.write("{headers}\n{contents}\n\ntime: {d_t}".format(headers=headers, contents=contents, d_t=d_t))


def main():
    lattice = generate_lattice()

    temperature_params = (MIN_TEMPERATURE, MAX_TEMPERATURE)
    cycle_params = (
        UPDATES_TO_SKIP,
        MAX_FLIP_ATTEMPTS,
        FLIPS_PER_UPDATE,
        UPDATES_PER_TEMPERATURE,
    )

    t_1 = datetime.datetime.now()

    results = simulate(temperature_params, cycle_params, J)

    t_2 = datetime.datetime.now()
    d_t = t_2 - t_1

    save_results(results, d_t)


main()
