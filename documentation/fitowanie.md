


# Fitowanie wyników symulacji

## Obróbka nagrań

Analizowane fragmenty nagrań to odcinki, w których zachodzi przejście fazowe.
Długości wszystkich fragmentów są takie same. Obejmują one czas przed przejściem fazowym,
w którym natężenie oklasków jest wysokie, oraz po przejściu fazowym,
kiedy nikt nie klaszcze - łączenie około 7 sekund.


Wszystkie nagrania zostały przekonwertowane na obrazy za pomocą programu *Sonic Visualizer*.
Pozwoliło to na wizualną analizę nagrań pod względem natężenia dźwięku, częstotliowości,
zakłóceń i innych - pozycja piksela w poziomie odpowiada momentowi na nagraniu,
pozycja w pionie odpowiada częstotliwości (w skali liniowej), a jasność piksela
to natężenie dźwięku. W celu odseparowania szumu, na obrazie rejestrowane były
tylko punkty o natężeniu powyżej pewnego progu.


*Sonic Visualiser* w akcji.

![Sonic Visualiser](https://github.com/micouy/ising/blob/master/documentation/sonic-visualiser.png)


Zakłócenie w nagraniu - dźwięki pianina nakładające
się na oklaski.

![Zakłócenie](https://github.com/micouy/ising/blob/master/documentation/artefact-1.png)


Jasności pikseli w każdej kolumnie obrazu zostały zsumowane, co odpowiadało
zsumowaniu natężenia każdej częstotliwości. Tak powstała lista zawierała tyle
wartości głośności, ile pikseli mieści się w szerokości obrazu. Została ona
podzielona na `n` równych części - odstępów w czasie. Wartość głośności w każdym
z tych kroków została uśredniona. Im większe `n`, tym więcej wartości zawiera lista
(ma większą "rozdzielczość"), ale również tym większe są odchylenia pomiędzy sąsiednimi wartościami,
co zakłóca obserwację ogólnej tendencji, np. spadku głośności.


Wartości średniej głośności w każdej części zostały znormalizowane, tj. podzielone
przez najwyższą z tych wartości. Taki zabieg ma kilka zalet:

* Nagrania są bardziej ujednolicone - wartości głośności mają wartość w przedziale
  od 0 do 1 w każdym nagraniu.
* Wartość magnetyzacji modelu Isinga również znajduje się przedziale od 0 do 1.
* Rzeczywista głośność nagrania nie jest ważna, ponieważ może różnić się zależnie od liczności
  publiczności i sprzętu, którym zostało wykonane. Badane jest przejście ze stanu, w którym
  klaszcze największa ilość osób (taki stan utrzymuje się przez względnie długi czas i wartość
  głośności nie ma dużych odchyleń) do stanu, w którym nikt nie klaszcze. Można
  nadać im wartości odpowiednio 1 oraz 0, a wszystkim pozostałym wartości pomiędzy.


Wykresy poniżej przedstawiają natężenie oklasków od czasu. Zostały wybrane
spośród wszystkich (dziewięciu) ze względu na największe podobieństwo do wykresów
magnetyzacji modelu Isinga dla `J = 1` i `k = 1`. Dla wykresów w tym dokumencie
`n = 40`, co odpowiada ilości wartości temperatury w symulacji.

![Chorzów 1](https://github.com/micouy/ising/blob/master/charts/chorzow-1.svg)
![Chorzów 5](https://github.com/micouy/ising/blob/master/charts/chorzow-5.svg)
![Chorzów 10](https://github.com/micouy/ising/blob/master/charts/chorzow-10.svg)
![Chorzów 12](https://github.com/micouy/ising/blob/master/charts/chorzow-12.svg)

Pozostałe wykresy natężenia oklasków można znaleźć [tutaj](https://docs.google.com/spreadsheets/d/1Xhf7dJrlBedGPSxnQW9vLhacw0Nkq-SSAyWvmPWIYgg/edit?usp=sharing).


## Symulacja modelu Isinga dla różnych wartości stałych J oraz k

Symulacja została przeprowadzona 16 razy - dla wszystkich par liczb `0.4`, `0.8`,
`1.2` oraz `1.6`. Wybrane przykłady:

![J = 1, k = 1](https://github.com/micouy/ising/blob/master/charts/simulation-J=1-k=1.svg)
![J = 0.8, k = 1.2](https://github.com/micouy/ising/blob/master/charts/simulation-J=0.8-k=1.2.svg)
![J = 0.8, k = 0.8](https://github.com/micouy/ising/blob/master/charts/simulation-J=0.8-k=0.8.svg)

Pozostałe wykresy magnetyzacji można znaleźć [tutaj](https://docs.google.com/spreadsheets/d/1doEWrnUNG1x7ro2eZrDH0M4yujCCvIWKe5AvqXeMVAI/edit?usp=sharing).
