#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>

struct Verbindung {

	Verbindung();
	Verbindung(size_t H, size_t C, size_t N, size_t O, size_t Si, size_t S, size_t Br, long double calculated_mass, long double deviation);

	// Die Variablen geben die Anzahl der Atome in der Verbindung an
	size_t H;
	size_t C;
	size_t N;
	size_t O;
	size_t Si;
	size_t S;
	size_t Br;

	long double calculated_mass;
	long double deviation;
	long double DBE;

	bool operator<(const Verbindung& b) const {
		return deviation < b.deviation;
	}

	bool operator>(const Verbindung& b) const {
		return deviation > b.deviation;
	}

	void print();

};

Verbindung::Verbindung()
	: H(0), C(0), N(0), O(0), Si(0), S(0), Br(0), calculated_mass(0), deviation(0), DBE(0)
{}

Verbindung::Verbindung(size_t _H, size_t _C, size_t _N, size_t _O, size_t _Si, size_t _S, size_t _Br, long double _calculated_mass, long double _deviation)
	: H(_H), C(_C), N(_N), O(_O), Si(_Si), S(_S), Br(_Br), calculated_mass(_calculated_mass), deviation(_deviation), DBE(0)
{
}

void Verbindung::print()
{
	std::cout << "Br = " << Br << "\tS = " << S << "\tSi = " << Si << "\tO = " << O << "\tN: " << N << "\tC = " << C << "\tH = " << H <<
		"\t Masse: " << calculated_mass << "; Abweichung: " << deviation << " DBE: " << DBE << std::endl;
}

long double double_bond_equivalent(Verbindung& A) {
	
	//DBE = [Summe über n_i(v_i - 2)]/2 + 1; n_i Anzahl eines Atoms; v_i Valenz des Atoms
	
	long double DBE = static_cast<long double>(A.H)	* (1 - 2)
					+ static_cast<long double>(A.C)	* (4 - 2)
					+ static_cast<long double>(A.N)	* (3 - 2)
					+ static_cast<long double>(A.O)	* (2 - 2)
					+ static_cast<long double>(A.Si)* (4 - 2)
					+ static_cast<long double>(A.S)	* (2 - 2)
					+ static_cast<long double>(A.Br)* (1 - 2);

	DBE = DBE / 2;
	DBE++;

	return DBE;
}

// Ist number eine ganze Zahl?
bool is_integer(long double number)
{
	// Teile number in eine ganze Zahl und einen Nachkommabereich
	// Siehe https://en.cppreference.com/cpp/numeric/math/modf

	long double vorkommastellen;
	long double nachkommastellen = std::modfl(number, &vorkommastellen);

	if (nachkommastellen == 0) {
		return true;
	}

	else {
		return false;
	}
}

// Ist number ein ganzzahliges vielfaches von 1/2 ?
bool is_half(long double number)
{
	// Teile number in eine ganze Zahl und einen Nachkommabereich
	// Siehe https://en.cppreference.com/cpp/numeric/math/modf

	long double vorkommastellen;
	long double nachkommastellen = std::modfl(number, &vorkommastellen);

	if (nachkommastellen == 0.5l) {
		return true;
	}

	else {
		return false;
	}
}

//void ergebniss_suche(std::vector<Verbindung> &alle_treffer, ) {
//
//}


int main(void) {

	// Masse, für die eine Summenformel gefunden werden soll
	long double observed_mass = 178.11;
	
	// Massen des durchschnittlichen Atoms der Elemente
	const long double mass_H	= (1.00784l + 1.00811l) / 2.0l;
	const long double mass_C	= (12.0096l + 12.0116l) / 2.0l;
	const long double mass_N	= (14.00643l + 14.00728l) / 2.0l;
	const long double mass_O	= (15.99903l + 15.99977l) / 2.0l;
	const long double mass_Si	= (28.084l + 28.086l) / 2.0l;
	const long double mass_S	= (32.059l + 32.076l) / 2.0l;
	const long double mass_Br	= (79.901l + 79.907l) / 2.0l;

	// Anzahl der sicher vorhandenen Atome
	const size_t Anzahl_H	= 0;
	const size_t Anzahl_C	= 0;
	const size_t Anzahl_N	= 0;
	const size_t Anzahl_O	= 0;
	const size_t Anzahl_Si	= 0;
	const size_t Anzahl_S	= 0;
	const size_t Anzahl_Br	= 0;


	long double tolerance = 0.0001l;// Maximal tolerierbare Abweichung

	size_t iteration_counter = 0;	// Zählt, wie viele Kombinationen getestet wurden
	size_t hits = 0;				// Zählt die Treffer mit einer Abweichung von höchstens tolerance

	// Grenzen der Iteration festlegen
	size_t H_max = static_cast<size_t>(ceil(observed_mass / mass_H));
	size_t C_max = static_cast<size_t>(ceil(observed_mass / mass_C));
	size_t N_max = static_cast<size_t>(ceil(observed_mass / mass_N));
	size_t O_max = static_cast<size_t>(ceil(observed_mass / mass_O));
	size_t Si_max = static_cast<size_t>(ceil(observed_mass / mass_Si));
	size_t S_max = static_cast<size_t>(ceil(observed_mass / mass_S));
	size_t Br_max = static_cast<size_t>(ceil(observed_mass / mass_Br));


	// Alternative manuelle Festlegung der Grenzen
	 //H_max = 60;
	 //C_max = 15;
	 //N_max = 5;
	 //O_max = 1;
	 //Si_max = 1;
	 //S_max = 5;
	 //Br_max = 5;

	std::cout << "Grenzen:" << std::endl
		<< "Br = " << Br_max << "\t"
		<< "S = " << S_max << "\t"
		<< "Si = " << Si_max << "\t"
		<< "O = " << O_max << "\t"
		<< "N = " << N_max << "\t"
		<< "C = " << C_max << "\t"
		<< "H = " << H_max << std::endl << std::endl;

	long double calculated_mass = 0.0l;	// Berechnete Masse der Summenformel
	long double deviation = 0.0;		// Abweichung zur gesuchten Summenformel
	long double DBE = 0;				// Double Bond equivalent

	std::vector<Verbindung> alle_treffer; // Alle zutreffenden Summenformeln werden hier gespeichert
	
	// Beginne Zeitmessung
	std::chrono::time_point Anfang = std::chrono::high_resolution_clock::now();

	for (size_t Br = Anzahl_Br; Br < Br_max; Br++) {
		for (size_t S = Anzahl_S; S < S_max; S++) {
			for (size_t Si = Anzahl_Si; Si < Si_max; Si++) {
				for (size_t O = Anzahl_O; O < O_max; O++) {
					for (size_t N = Anzahl_N; N < N_max; N++) {
						for (size_t C = Anzahl_C; C < C_max; C++) {
							for (size_t H = Anzahl_H; H < H_max; H++) {
							
								iteration_counter++;

								// berechne Masse
								calculated_mass = static_cast<long double>(Br) * mass_Br
												+ static_cast<long double>(S) * mass_S
												+ static_cast<long double>(Si) * mass_Si
												+ static_cast<long double>(O) * mass_O
												+ static_cast<long double>(N) * mass_N
												+ static_cast<long double>(C) * mass_C
												+ static_cast<long double>(H) * mass_H;

								// berechne Abweichung
								deviation = abs(observed_mass - calculated_mass) / observed_mass;

								// vergleiche Abweichung und Toleranz
								if (deviation < tolerance) {

									Verbindung A(H, C, N, O, Si, S, Br, calculated_mass, deviation);

									// berechne Doppelbindungsäquivalente
									DBE = double_bond_equivalent(A);

									// Entscheide, ob die Summenformel aufgrund der Doppelbindungsäquivalente sinnvoll ist
									if (is_integer(DBE) && DBE >= 0) {
										A.DBE = DBE;
										alle_treffer.emplace_back(A);
										hits++;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	std::chrono::time_point Ende = std::chrono::high_resolution_clock::now();
	auto Rechenzeit = std::chrono::duration_cast<std::chrono::milliseconds>(Ende - Anfang);

	std::vector<Verbindung> alle_treffer_sorted = alle_treffer;
	std::sort(alle_treffer_sorted.begin(), alle_treffer_sorted.end());

	for (auto& x : alle_treffer_sorted) {
		x.print();
	}

	std::cout  << "Es wurden " << static_cast<long double>(iteration_counter) / (1000.0l * 1000.0l * 1000.0l) << " Milliarden Kombinationen in " << Rechenzeit.count() << " Millisekunden getestet; Treffer: " << hits << std::endl;
	std::cout  << static_cast<size_t>((static_cast<long double>(iteration_counter) / (static_cast<long double>(Rechenzeit.count() * 1000)))) << " Millionen berechnungen pro Sekunde" << std::endl;;
	return 0;
}


	