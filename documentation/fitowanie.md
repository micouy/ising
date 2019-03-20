


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


Wykresy poniżej przedstawiają natężenie oklasków od czasu na nagraniach oraz uśrednioną i znormalizowaną
głośność oklasków ze wszystkich nagrań.

![Oklaski](https://github.com/micouy/ising/blob/master/charts/recordings/recordings.svg)


## Symulacja modelu Isinga dla różnych wartości stałych `J` oraz `h`

Wykresy poniżej przedstawiają wartość magnetyzacji w danej temperaturze. Wartości `h` serii
są podane w legendach wykresów.

![Symulacja](https://github.com/micouy/ising/blob/master/charts/simulation/results.svg)


Pozostałe wyniki można zobaczyć [tutaj](https://docs.google.com/spreadsheets/d/1doEWrnUNG1x7ro2eZrDH0M4yujCCvIWKe5AvqXeMVAI/edit?usp=sharing).
