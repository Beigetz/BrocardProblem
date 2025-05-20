#include <iostream>
#include <mpirxx.h> // or #include <gmpxx.h> MPIR Library's header or GMP Library's header

mpz_class multifactorial(mpz_class n, mpz_class m);
mpz_class calculate_b(int p, int j, int m, int n);
mpz_class calculate_delta(int p, int j, int m, int n);
mpz_class calculate_a(int p, int j, int m, int n);
mpz_class calculate_c(int p, int j, int m, int n);

// GCD of 2 numbers
mpz_class gcd(mpz_class x, mpz_class y) {
	x = abs(x);
	y = abs(y);
	if (y == 0) return x;
	mpz_class result;
	mpz_gcd(result.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
	return result;
}

// GCD of 3 numbers 
mpz_class gcd3(mpz_class a, mpz_class b, mpz_class c) {
	return gcd(gcd(a, b), c);
}

// Verify if the ratio of square roots is coprime
void check_square_ratio(mpz_class x, mpz_class y, mpz_class z) {
	mpz_class a = sqrt(x);
	mpz_class b = sqrt(y);
	mpz_class c = sqrt(z);

	// original integer ratio
	mpz_class g = gcd3(x, y, z);
	mpz_class rx = x / g;
	mpz_class ry = y / g;
	mpz_class rz = z / g;

	std::cout << "Original numbers: " << x << ", " << y << ", " << z << "\n";
	std::cout << "Ratio: " << rx << ":" << ry << ":" << rz << "\n";


	// The ratio of square roots
	mpz_class g_sqrt = gcd3(a, b, c);
	mpz_class ra = a / g_sqrt;
	mpz_class rb = b / g_sqrt;
	mpz_class rc = c / g_sqrt;

	std::cout << "Square roots: " << a << ", " << b << ", " << c << "\n";
	std::cout << "Ratio of square roots: " << ra << ":" << rb << ":" << rc << "\n";

	// Coprime or not
	if (g_sqrt == 1) {
		std::cout << "The ratio of square roots is coprime." << "\n";
	}
	else {
		std::cout << "The ratio of square roots is NOT coprime (GCD = " << g_sqrt << ")." << "\n";
	}
}

int main() {
	int p, j, m, n, q;
	q = 1;
	while (q) {
		std::cout << "input parameters p, j,m, and n, exclude (j,m,n)=(2,1,-1):\n";
		std::cout << "input positive integer p>=1: ";
		std::cin >> p;
		std::cout << "input positive integer j>=2: ";
		std::cin >> j;
		std::cout << "input positive integer m>=1: ";
		std::cin >> m;
		std::cout << "input integer n>=-1: ";
		std::cin >> n;

		mpz_class term_2jm_plus_n_mpz = (mpz_class(m) << j) + n;
		mpz_class c_val_mpz = calculate_c(p, j, m, n);
		mpz_class c_val_abs = abs(c_val_mpz);
		mpz_class a_val_mpz = calculate_a(p, j, m, n);
		mpz_class multi = multifactorial((m << j) + n, m);
		mpz_class b_val_mpz = calculate_b(p, j, m, n);

		// Calculate left-hand side
		mpz_class lhs1 = multi;
		mpz_class lhs2 = multi;
		for (int i = 1; i < p; ++i) { //multiplying `lhs2` by `multi for (p-1) times
			lhs2 *= multi;
		}
		mpz_class lhs = lhs2 + c_val_abs * c_val_abs; // left-hand side

		// Calculate right-hand side
		mpz_class a_minus_c_mpz = a_val_mpz - c_val_mpz;
		mpz_class rhs = a_minus_c_mpz * a_minus_c_mpz; // right-hand side
		// mpz_class del = calculate_delta(p, j, m, n);

		std::cout << "p: " << p << ", ";
		std::cout << "j: " << j << ", ";
		std::cout << "m: " << m << ", ";
		std::cout << "n: " << n << "\n";
		/*
		std::cout << "a:" << a_val_mpz << ",\n b:" << b_val_mpz << ",\n c:" << c_val_abs << ",\n (2^j*m+n)!_m" << multi << ",\n delta:" << del << "\n";
		std::cin >> q;
		*/
		if ((j == 2 && m == 1 && n == -1) || (p < 1) || (j < 2) || (m < 1) || (n < -1)) {
			std::cout << "Error\n";
		} else if (lhs==rhs) {
		
			if (p == 1) {
				if (m==1) {
					if (j == 2) {
						if (n == 0 || n == 1) {
							std::cout << "The Brocard number : (" << (m<<j)+n << ", " << sqrt(rhs) << ")\n" ;
							std::cout << "(" << term_2jm_plus_n_mpz << ")!  + " << c_val_abs <<  " = " << a_minus_c_mpz << " ^ 2\n";
						}
					} else std::cout << "(" << term_2jm_plus_n_mpz << ")!  + " << c_val_abs << " ^ 2" << " = " << a_minus_c_mpz << " ^ 2\n";
					else if (j == 3) {
						if (n==-1) {
							std::cout << "The Brocard number : (" << (m << j) + n << ", " << sqrt(rhs) << ")\n";
							std::cout << "(" << term_2jm_plus_n_mpz << ")!  + " << c_val_abs <<  " = " << a_minus_c_mpz << " ^ 2\n";
						}
					}else std::cout << "(" << term_2jm_plus_n_mpz << ")!  + " << c_val_abs << " ^ 2" << " = " << a_minus_c_mpz << " ^ 2\n";
				}
				else {
					std::cout << "(" << term_2jm_plus_n_mpz << ")!_" << m  << "  + (" << c_val_abs << ")^ 2" << " = (" << a_minus_c_mpz << ")^2\n";
				}
			}
			else {
				if (p == 2) {
					/*mpz_class a = multifactorial(term_2jm_plus_n_mpz, m);
					mpz_class b = c_val_abs;
					mpz_class c = a_minus_c_mpz;*/
					check_square_ratio(multifactorial(term_2jm_plus_n_mpz, m), c_val_abs, a_minus_c_mpz);
				}
				if (m == 1) {
					std::cout << "(" << term_2jm_plus_n_mpz << ")!^" << p  << "  + (" << c_val_abs << ")^ 2" << " = (" << a_minus_c_mpz << ")^2\n";
				}
				else {
					std::cout << "((" << term_2jm_plus_n_mpz << ")!_" << m << "^" << p << " + (" << c_val_abs << ")^2" << " = (" << a_minus_c_mpz << ")^2\n";
				}
			}
			// std::cout << term_2jm_plus_n_mpz << "\n";
			std::cout << "Left hand side value : " << lhs << "\n";
			std::cout << "Right hand side value : " << rhs << "\n";
			std::cout << "Square root of value : " << sqrt(rhs) << "\n";
			std::cout << "difference : " << lhs - rhs << "\n";
		}
		else {
			break;
		}
		/*mpz_clear(term_2jm_plus_n_mpz);*/
		std::cout << "\n continue? yes:1, no:0 :";
		std::cin >> q;
	}
	return 0;
}
// multifactrial function
mpz_class multifactorial(mpz_class n, mpz_class m) {
	if (n <= 0) return 1;  // if n<=0 1
	mpz_class result = 1;

	// Calculate the product while decreasing from n
	for (mpz_class i = n; i > 0; i -= m) {
		result *= i;
	}

	return result;
}
// B(p, j, m, n) function
mpz_class calculate_b(int p, int j, int m, int n) {
	mpz_class product = 1;
	int limit = 1 << (j - 2);
	for (int k = 0; k < limit; ++k) {
		mpz_class term1 = ((mpz_class(1) << (j - 1)) + 2 * k + 1) * m + n;
		mpz_class term2 = (mpz_class((2 * k + 2) * m + n));
		product *= term1 * term2;
	}
	mpz_class result = product;
	if (p > 1) {
		for (int i = 1; i < p; ++i) {
			result *= product;
		}
	}
	return result;
}
// Delta(p,j,m,n) function
mpz_class calculate_delta(int p, int j, int m, int n) {
	if (m != 1 || n != -1) return 0;
	mpz_class d = multifactorial((1 << j) + n, 1);

	if (p > 1) {
		for (int i = 1; i < p; ++i) {
			d *= d;
		}
	}
	mpz_class delta = d * ((m == 1) ? 1 : 0)*((n == -1) ? 1 : 0) / calculate_b(p, j, m, n);
	//mpz_class result = delta;

	return delta;
}
// A(p, j, m, n) function
mpz_class calculate_a(int p, int j, int m, int n) {
	mpz_class product = 1;
	int limit = 1 << (j - 2); // 2^(j-2)

	// Calculate m-th multifactorial of n
	mpz_class n_multifac_mpz = multifactorial(n <= 0 ? 1 : n, m);

	// Prod_{k=0}^{2^{j-2}-1}
	for (int k = 0; k < limit; ++k) {
		mpz_class term1 = ((mpz_class(1) << (j - 1)) + 2 * k + 2) * m + n;
		mpz_class term2 = ((mpz_class(k) << 1) + 1) * m + n;
		product *= term1 * term2;
	}

	// p-th power
	mpz_class term = product * n_multifac_mpz;
	mpz_class result = 1;

	for (int i = 0; i < p; ++i) {
		result *= term;
	}

	result += calculate_delta(p, j, m, n); // adding delta

	return result;
}
// C(p, j, m, n) function
mpz_class calculate_c(int p, int j, int m, int n) {
	return (calculate_a(p, j, m, n) - calculate_b(p, j, m, n)) / 2;
}
