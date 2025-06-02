#pragma once
#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <functional>
#include <fstream>

// Sta�e globalne
extern const double P = 1000.0;                    // [N] Si�a dzia�aj�ca
extern const double E = 207e9;                     // [Pa] Modu� Younga
extern const double rho = 7800.0;                  // [kg/m^3] G�sto�� materia�u
extern const double u_max = 5e-3;                  // [m] Maksymalne ugi�cie
extern const double sigma_max = 300e6;             // [Pa] Maksymalne napr�enie
extern const double PI = 3.14159265358979323846;   // Liczba PI

extern int function_call_count = 0;                // Licznik wywo�a� funkcji celu
extern const double l_min_g = 0.2;  // [m]
extern const double l_max_g = 1.0;  // [m]
extern const double d_min_g = 0.01; // [m]
extern const double d_max_g = 0.05; // [m]

// ============================   [ Metody obliczeniowe ]   ============================
// Znajdowanie przedzia�u wyst�powania minimum funkcji liniowej - metoda ekspansji
std::pair<double, double> metoda_ekspansji(const std::function<double(double)>& f, double x1, double x2, double alfa, int N_max);

// Znajdowanie minimum funkcji liniowej - metoda z�otego podzia�u
double metoda_zlotego_podzialu(const std::function<double(double)>& f, double a, double b, int N_max, const double epsilon);
// Znajdowanie optymalnego minimum funkcji - metod� Powella
std::vector<double> metoda_Powella_v2(const std::function<double(double, std::vector<double>, double)>& f, std::vector<double> x_0, double w, double c, int N_max, const double epsilon);
// =====================================================================================

// ========================   [ Funkcje problemu testowego ]   =========================

// Funkcja pierwsza, wg dostarczonego wzoru
double f1(double x1, double x2, double a);

// Funkcja druga, wg dostarczonego wzoru
double f2(double x1, double x2, double a);

// Funkcja wa�ona problemu testowego
double problem_jednokryterialny2(double w, std::vector<double> x, double a);

// Wywo�anie optymalizacji dla problemu rzeczywistego
void problem_testowy();

// =====================================================================================

// ======================   [ Funkcje problemu rzeczywistego ]   =======================

// Funkcja obliczaj�ca mas�
double f_m(double l, double d);

// Funkcja obliczaj�ca ugi�cie
double f_u(double l, double d);

// Funkcja obliczaj�ca napr�enie
double sigma(double l, double d);

// Sprawdzanie czy s� spe�nione zadane ograniczenia
bool sprawdzanie_ograniczen(double l, double d);

// Funkcja obliczaj�ca kar� zewn�trzn�
double funkcja_kary(double l, double d, double c);

// Funkcja wa�ona problemu rzeczywistego
double funkcja_wazona(double w, std::vector<double> x, double c);

// Wywo�anie optymalizacji dla problemu rzeczywistego
void problem_rzeczywisty();

// =====================================================================================

// ========================   [ Dodatkowe metody pomocnicze ]   ========================

// Generowanie wektora o d�ugo�ci ilosc ze zmiennoprzecinkowymi warto�ciami losowymi z przedzia�u [a, b]
std::vector<double> generowanie_wektora_wartosci_losowych(double a, double b, int ilosc);

// Generowanie zmiennoprzecinkowej warto�ci losowej z przedzia�u [a, b]
double generowanie_wartosci_losowej(double a, double b);

// =====================================================================================