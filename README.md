# Optymalizacja Wielokryterialna – Projekt C++

## 📌 Cel projektu

Celem projektu jest zapoznanie się z metodami optymalizacji wielokryterialnej oraz wyznaczenie rozwiązań minimalnych w sensie Pareto. Projekt obejmuje zarówno testową funkcję celu, jak i problem rzeczywisty (optymalizacja kształtu belki).

## ⚙️ Zastosowane metody

- **Metoda kryterium ważonego** – przekształcenie problemu wielokryterialnego na jednokryterialny:
  $f(\mathbf{x}) = w \cdot f_1(\mathbf{x}) + (1 - w) \cdot f_2(\mathbf{x}), \quad w \in [0, 1]$
- **Metoda Powella** – do znajdowania minimum funkcji celu.
- **Złoty podział** – do znajdowania minimum na kierunku.
- **Metoda ekspansji** – do wyznaczania początkowego przedziału.
- **Zewnętrzna funkcja kary** – do uwzględnienia ograniczeń w problemie rzeczywistym.

## 🧪 Część 1 – Testowa funkcja celu

Zdefiniowano dwie funkcje celu:
- $f_1(x_1, x_2) = a \cdot ((x_1 - 2)^2 + (x_2 - 2)^2)$
- $f_2(x_1, x_2) = \frac{1}{a} \cdot ((x_1 + 2)^2 + (x_2 + 2)^2)$

Gdzie:
- $a \in \{1, 10, 100\}$
- $x_1^{(0)}, x_2^{(0)} \in [-10, 10]$ 

Dla każdej wartości parametru $a$, wykonano **101 optymalizacji** dla $w = \{0, 0.01, ..., 1\}$ Wyniki zapisano w formacie `.xlsx`, a rozwiązania Pareto przedstawiono na wykresach

## 🏗️ Część 2 – Problem rzeczywisty (belka)

Optymalizujemy belkę o przekroju kołowym poddaną obciążeniu $P = 1\ \text{kN}$. Zmienne decyzyjne:
- długość belki $l \in [200, 1000]\ \text{ mm}$
- średnica $d \in [10, 50]\ \text{ mm}$

### Funkcje celu:
- **Masa belki (f₁)**:\
  $f_1(l, d) = \rho \cdot \pi \cdot \left( \frac{d}{2} \right)^2 \cdot l$
- **Ugięcie belki (f₂)**:\
  $f_2(l, d) = \frac{64 \cdot P \cdot l^3}{3 \cdot E \cdot \pi \cdot d^4}$

### Ograniczenia:
- Maksymalne ugięcie $u \leq 5\ \text{ mm}$
- Maksymalne naprężenie $\sigma \leq 300\ \text{ MPa}$

### Naprężenie:
$\sigma = \frac{32 \cdot P \cdot l}{\pi \cdot d^3}$

Wykonano 101 optymalizacji dla $w \in [0, 1]$. Wyniki zapisano w pliku `.xlsx`, a rozwiązania Pareto przedstawiono na wykresach.

## 📊 Wyniki

Wyniki zostały przedstawione w formie:
- tabel `.xlsx`
- wykresów Pareto dla każdej wartości parametru `a` oraz dla problemu rzeczywistego

## 📄 Sprawozdanie

Dokumentacja w formacie `.docx`/`.pdf` zawiera:
- opis metod
- analizę wyników
- wnioski
- kod źródłowy
- wyniki w `.xlsx`

## 👨‍💻 Autor

Projekt wykonany w ramach zajęć z optymalizacji, AGH WIMiIP  
Autor: Karol Woda
