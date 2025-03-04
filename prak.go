package main

import (
	"bufio"
	"errors"
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"os"
	"reflect"
	"strconv"
	"strings"
	"time"
)

// gcd calculates the greatest common divisor of two integers a and b
// using the Euclidean algorithm. It returns the largest integer that
// divides both a and b without leaving a remainder.

func gcd(a, b int64) int64 {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// lcm calculates the least common multiple of two integers a and b.
// It returns the smallest positive integer that is a multiple of both a and b.
func lcm(a, b *big.Int) *big.Int {
	if a.Cmp(big.NewInt(0)) == 0 || b.Cmp(big.NewInt(0)) == 0 {
		return big.NewInt(0)
	}
	gcdVal := new(big.Int).GCD(nil, nil, a, b)
	result := new(big.Int).Mul(new(big.Int).Div(a, gcdVal), b)
	return result
}

// Power calculates a to the power b modulo m.
// It uses the big package to compute the result.
// The arguments a, b and m are int64 values.
func Power(a, b, m int64) int64 {
	bigA := big.NewInt(a)
	bigB := big.NewInt(b)
	bigM := big.NewInt(m)
	result := new(big.Int).Exp(bigA, bigB, bigM)
	return result.Int64()
}

// legendreSymbol returns the Legendre symbol of a with respect to p.
// It returns 0 if a is 0, 1 if a is a quadratic residue modulo p,
// and -1 if a is a quadratic nonresidue modulo p.
func legendreSymbol(a, p int64) int64 {
	a %= p
	if a == 0 {
		return 0
	}
	if Power(a, (p-1)/2, p) == 1 {
		return 1
	}
	return -1
}

// factorize returns the prime factorization of n in the finite field of integers modulo p.
// The returned slice contains the prime factors of n in the order in which they were found.
func factorize(n int64, p int64) []int64 {
	if n <= 1 {
		return []int64{} // Обработка чисел меньше или равных 1
	}

	factors := []int64{}

	// Обработка делимости на 2
	for ; n%2 == 0; n /= 2 {
		factors = append(factors, 2)
	}

	// Проверка делимости на нечётных числах до sqrt(n)
	var i_help big.Int
	for i := int64(3); i_help.Exp(big.NewInt(i), big.NewInt(2), big.NewInt(p)).Int64() <= n; i += 2 { //i*i <= n; i += 2 {
		for ; n%i == 0; n /= i {
			factors = append(factors, i)
		}
	}

	// Если n > 1 после цикла, значит n - простое число
	if n > 1 {
		factors = append(factors, n)
	}

	return factors
}

// In checks if target is in slice.
// It returns true if target is in slice, and false otherwise.
// The function uses the DeepEqual function from the reflect package to compare the elements of slice with target.
func In(slice [][]int64, target []int64) bool {
	for i := 0; i < len(slice); i++ {
		if reflect.DeepEqual(slice[i], target) {
			return true
		}
	}
	return false
}

// In_2 checks if target is in slice.
// It returns true if target is in slice, and false otherwise.
// The function uses the DeepEqual function from the reflect package to compare the elements of slice with target.
func In_2(slice []int64, target int64) bool {
	for i := 0; i < len(slice); i++ {
		if reflect.DeepEqual(slice[i], target) {
			return true
		}
	}
	return false
}

// Index returns the index of target in slice, or -1 if target is not found in slice.
// The function uses the DeepEqual function from the reflect package to compare the elements of slice with target.
func Index(slice [][]int64, target []int64) int64 {
	var index int64
	var Length int64
	Length = int64(len(slice))
	for index = 0; index < Length; index++ {
		if reflect.DeepEqual(slice[index], target) {
			break
		}
	}
	return index
}

type simple_field struct {
	p              int64
	field_elements []int64
	squares        map[int64]int64
}

// New_Field returns a pointer to a simple_field with the value of p being 0,
// field_elements being [0], and squares being an empty map.
// The map squares does not contain the key 0.
func New_Field() *simple_field {
	return &simple_field{
		p:              0,
		field_elements: []int64{0},
		squares:        map[int64]int64{}, //убрал 0: 0
	}
}

// Symmetrical_field_construction constructs a simple field symmetrically by
// filling in the field_elements slice with values from 1 to p/2 (inclusive)
// and populating the squares map with the squares of the values in
// field_elements, with the index as the key and the value as the value.
func (field *simple_field) Symmetrical_field_construction() {
	var limit int64 = int64(field.p/2 + 1)
	for i := int64(1); i < limit; i++ {
		field.field_elements = append(field.field_elements[:], i)
	}
	for index, value := range field.field_elements {
		square := value*value%field.p + field.p%field.p
		field.squares[square] = int64(index)
	}
}

// Asymmetrical_field_construction constructs a simple field asymmetrically by
// filling in the field_elements slice with values from 1 to p-1 (inclusive).
func (field *simple_field) Asymmetrical_field_construction() {
	for i := int64(1); i < field.p; i++ {
		field.field_elements = append(field.field_elements[:], i)
	}
}

// The_opposite_in_symmetric_field finds the opposite of the element 'a' in a
// symmetric field using the extended Euclidean algorithm. It first reduces
// 'a' modulo 'p' to ensure that it is an element of the symmetric field.
// The algorithm works by setting up a system of congruences, where each
// congruence is of the form 'a' * x + 'b' * y = 'n' (mod 'p'). The algorithm
// then successively replaces 'a' and 'b' with 'b' and 'r' (the remainder of
// 'a' divided by 'b'), and 'n' with 'r', until 'r' becomes 0. The value of
// 'y' is then the opposite of 'a' in the symmetric field.
func (field *elliptical_curve) The_opposite_in_symmetric_field(a int64) int64 {
	a = (a%field.p + field.p) % field.p // Точно элемент ассиметричного поля
	var r, y1, y2 int64
	r, y1, y2 = 1, 1, 0
	n := field.p
	n_old := n
	for i := 1; i > 0; i++ {
		if r != 0 {
			q := n / a
			r = n - q*a
			y := y2 - q*y1
			n = a
			a = r
			y2 = y1
			y1 = y
		} else {
			break
		}
	}
	return (y2%n_old + n_old) % n_old
}

type elliptical_curve struct {
	a              int64
	b              int64
	p              int64
	field_elements [][]int64
}

// Point_creator calculates the points on the elliptic curve of the form
// y^2 = x^3 + ax + b (mod p) for a given x-coordinate. It returns a slice
// of two points, each of which is a slice of two int64 values. If the point
// is an element of the asymmetric field, it will return the point twice.
// If the point is not an element of the asymmetric field, it will return an
// empty slice.
func (field *elliptical_curve) Point_creator(x int64) [][]int64 {
	var Points [][]int64
	var Y_square int64
	Y_square = ((field.power(x, 3)+field.a*x+field.b)%field.p + field.p) % field.p // Поменял
	if Y_square == 0 {
		Points = append(Points, []int64{x, 0})
		Points = append(Points, []int64{x, 0})
		return Points
	}
	if reflect.DeepEqual(field.Square_root(Y_square), []int64{-1, -1}) {
		return [][]int64{{0}}
	} else {
		Points = append(Points, []int64{x, field.Square_root(Y_square)[0]})
		Points = append(Points, []int64{x, field.Square_root(Y_square)[1]})
	}
	return Points
}

// Square_root finds the square root of a given number in the field of the elliptic curve. If the number does not have a square root, it returns [-1, -1].
func (field *elliptical_curve) Square_root(k int64) []int64 {
	var b int64
	var s int64
	if k == 0 {
		return []int64{0, 0}
	}
	if legendreSymbol(k, field.p) == -1 {
		return []int64{-1, -1}
	} else {
		if legendreSymbol(-1, field.p) == -1 {
			b = -1
		} else {
			b = 2
			for legendreSymbol(b, field.p) != -1 {
				b += 1
			}
		}
		s = 0
		p_1 := field.p - 1
		for p_1%2 == 0 {
			p_1 = p_1 / 2
			s += 1
		}
		var t_help big.Int
		t_help.Exp(big.NewInt(2), big.NewInt(s), big.NewInt(field.p))
		t_help.Div(big.NewInt(field.p-1), &t_help)
		t := t_help.Int64()
		reversed_k := ((field.Asymmetry_to_symmetry(field.The_opposite_in_symmetric_field(k)) % field.p) + field.p) % field.p
		C_0 := (field.power(b, t)%field.p + field.p) % field.p
		r := (field.power(k, (t+1)/2)%field.p + field.p) % field.p
		for i := int64(1); i < s; i++ {
			var d_i int64
			bigP := big.NewInt(field.p)
			bigR := big.NewInt(r)
			bigReversedK := big.NewInt(reversed_k)
			bigMul := new(big.Int).Mul(bigR, bigR)
			bigMul.Mul(bigMul, bigReversedK)
			bigExp := new(big.Int).Exp(big.NewInt(2), big.NewInt(s-i-1), nil)
			bigDI := new(big.Int).Exp(bigMul, bigExp, bigP)
			bigDI.Mod(bigDI.Add(bigDI, bigP), bigP)
			if bigDI.IsInt64() {
				d_i = bigDI.Int64()
			} else {
				panic("Overflow detected in d_i computation")
			}
			if d_i%field.p == field.p-1 {
				var r_help big.Int
				r_help.Mul(big.NewInt(r), big.NewInt(C_0))
				r_help.Mod(&r_help, big.NewInt(field.p))
				r = r_help.Int64()
			}
			C_0 = (field.power(C_0, 2)%field.p + field.p) % field.p
		}
		return []int64{field.Asymmetry_to_symmetry(r), -(field.Asymmetry_to_symmetry(r))}
	}
}

// elliptical_curve_construction constructs the elliptic curve of the form y^2 = x^3 + ax + b (mod p)
// for a given simple field. It returns the length of the resulting field. The function first checks that
// the smoothness condition is met, and then calculates the points on the elliptic curve. For each point,
// it adds the point to the field, and then adds its opposite in the symmetric field to the field if the
// point is not the point at infinity.
func (field *elliptical_curve) elliptical_curve_construction(Simple_field simple_field) int64 {
	var Y_square int64
	field.p = Simple_field.p
	if ((-4*int64(math.Pow(float64(field.a), float64(3)))-27*field.b*field.b)%Simple_field.p+Simple_field.p)%Simple_field.p == 0 {
		fmt.Println("The smoothness condition is not met")
		os.Exit(1)
	}
	for i := int64(0); i < Simple_field.p; i++ {
		x := float64(i)
		Y_square = ((int64(math.Pow(x, float64(3)))+field.a*int64(x)+field.b)%Simple_field.p + Simple_field.p) % Simple_field.p
		value, ok := Simple_field.squares[Y_square]
		if x > float64(Simple_field.p/2) {
			x = x - float64(Simple_field.p)
		}
		element := []int64{int64(x), value}
		if ok {
			if Y_square == 0 {
				field.field_elements = append(field.field_elements, element)
			} else {
				field.field_elements = append(field.field_elements, element)
				element := []int64{int64(x), -value}
				field.field_elements = append(field.field_elements, element)
			}
		} else {
			continue
		}
	}
	return int64(len(field.field_elements))
}

// power calculates a^k in the field of the elliptic curve.
// It takes advantage of the exponentiation by squaring algorithm to
// reduce the number of multiplications required.
func (field *elliptical_curve) power(a int64, k int64) int64 {
	K_binary := strconv.FormatInt(int64(k), 2)
	Rune_k := []rune(K_binary)
	var b int64
	b = 1
	a = (a%field.p + field.p) % field.p
	var A_help big.Int
	var b_help big.Int
	if k == 0 {
		return b
	}
	A := a
	if string(Rune_k[len(K_binary)-1]) == "1" {
		b = a
	}
	for i := len(K_binary) - 2; i >= 0; i-- {
		A_help.Exp(big.NewInt(A), big.NewInt(2), big.NewInt(field.p))
		A = (A_help.Int64()%field.p + field.p) % field.p
		if string(Rune_k[i]) == "1" {
			b_help.Mul(big.NewInt(A), big.NewInt(b))
			b_help.Mod(&b_help, big.NewInt(field.p))
			b = b_help.Int64()
		}
	}
	b = b % field.p
	return b

}

// Adding_points adds two points on the elliptic curve. If the points are the same,
// it doubles the point. If the points are different, it adds them. If the points are
// the same and the y-coordinate of the point is 0, it returns the point at infinity.
// If the points are different but have the same x-coordinate, it returns the point
// at infinity. Otherwise, it returns the sum of the two points. The function takes
// advantage of the property of the elliptic curve that the sum of three points is
// the same as the sum of two points, and the property that the sum of a point and
// its opposite is the point at infinity.
func (field *elliptical_curve) Adding_points(Points [][]int64) [][]int64 {
	var X3 int64 = 0
	var Y3 int64 = 0
	var Result [][]int64 = [][]int64{{X3, Y3}}
	var x_3_help, y3_help big.Int
	if len(Points) == 1 || reflect.DeepEqual(Points[0], Points[1]) {
		if reflect.DeepEqual(Points[0], []int64{0}) || Points[0][1] == 0 {
			return [][]int64{{0}}
		} else {
			rev := field.The_opposite_in_symmetric_field(2 * Points[0][1])
			x_3_help.Mul(big.NewInt(3*field.power(Points[0][0], 2)+field.a), big.NewInt(rev))
			x_3_help.Mod(&x_3_help, big.NewInt(field.p))
			X_3 := x_3_help.Int64()
			X3 = field.Asymmetry_to_symmetry(x_3_help.Int64()*x_3_help.Int64() - 2*Points[0][0])
			y3_help.Mul(big.NewInt(X_3), big.NewInt(Points[0][0]-X3))
			Y3 = field.Asymmetry_to_symmetry(y3_help.Int64() - Points[0][1])
		}
	} else {
		if reflect.DeepEqual(Points[0], []int64{0}) {
			Result = [][]int64{Points[1]}
			return Result
		}
		if reflect.DeepEqual(Points[1], []int64{0}) {
			Result = [][]int64{Points[0]}
			return Result
		}
		if Points[0][0] == Points[1][0] {
			return [][]int64{{0}}
		}
		rev := field.The_opposite_in_symmetric_field(Points[1][0] - Points[0][0])
		x_3_help.Mul(big.NewInt(Points[1][1]-Points[0][1]), big.NewInt(rev))
		x_3_help.Mod(&x_3_help, big.NewInt(field.p))
		X_3 := x_3_help.Int64()
		X3 = field.Asymmetry_to_symmetry(x_3_help.Int64()*x_3_help.Int64() - Points[0][0] - Points[1][0])
		y3_help.Mul(big.NewInt(X_3), big.NewInt(Points[0][0]-X3))
		Y3 = field.Asymmetry_to_symmetry(y3_help.Int64() - Points[0][1])
	}
	Result = [][]int64{{field.Asymmetry_to_symmetry(X3), field.Asymmetry_to_symmetry(Y3)}}
	return Result
}

// New_elliptical_curve creates a new elliptic curve with default values for a, b, p, and field_elements.
// The default values are a = 0, b = 0, p = 0, and field_elements = {{0}}.
// The function returns a pointer to the newly created elliptic curve.
func New_elliptical_curve() *elliptical_curve {
	return &elliptical_curve{
		a:              0,
		b:              0,
		p:              0,
		field_elements: [][]int64{{0}},
	}
}

// Asymmetry_to_symmetry takes an element a of the asymmetric field and returns its
// representation in the symmetric field. The function first reduces 'a' modulo 'p'
// to ensure that it is an element of the symmetric field. It then subtracts 'p' from
// 'a' if 'a' is greater than (p-1)/2. The function returns the result.
func (field *elliptical_curve) Asymmetry_to_symmetry(a int64) int64 {
	a = (a%field.p + field.p) % field.p
	if a > (field.p-1)/2 {
		a = a - field.p
	}
	return a
}

// Calculating_point_of_given_multiplicity calculates the point P*Q on the elliptic curve where
// Q is the multiplicity and P is the point given as argument. It first reduces the multiplicity
// modulo the order of the elliptic curve. It then uses the double and add algorithm to
// calculate the point P*Q. The function returns the result as a slice of two int64 values.
func (field *elliptical_curve) Calculating_point_of_given_multiplicity(Point [][]int64, multiplicity int64) [][]int64 {
	Result := [][]int64{{0}}
	pow := multiplicity
	Tmp := [][]int64{{Point[0][0], Point[0][1]}}
	if pow < 0 {
		pow = -pow
		Tmp[0][1] = -Tmp[0][1]
	}
	if pow == 0 || reflect.DeepEqual(Tmp, [][]int64{{0}}) {
		return [][]int64{{0}}
	}
	for pow != 0 {
		if pow%2 == 1 {
			Result = field.Adding_points([][]int64{Result[0], Tmp[0]})
		}
		pow = pow / 2
		Tmp = field.Adding_points(Tmp)
	}
	return Result
}

// Big_step_for_giant_small_step_for_baby calculates the order of the point P using the big step/giant step algorithm.
// The function takes a point on the elliptic curve as argument and returns the order of the point as an int64 value.
func (field *elliptical_curve) Big_step_for_giant_small_step_for_baby(Point [][]int64) int64 {
	Q := field.Calculating_point_of_given_multiplicity([][]int64{Point[0]}, field.p+1)
	m := int64(math.Pow(float64(field.p), 0.25)) + 1
	var jPs [][]int64
	var limit int64 = m + 1
	for j := int64(0); j < limit; j++ {
		jPs = append(jPs, field.Calculating_point_of_given_multiplicity([][]int64{Point[0]}, j)[0])
	}
	var l int64
	var k int64
	for k = -m; k < limit; k++ {
		help_degree := new(big.Int).Mul(big.NewInt(2), big.NewInt(m))
		help_degree.Mul(help_degree, big.NewInt(k))
		degree := help_degree.Int64()
		Q_ := field.Adding_points([][]int64{Q[0], field.Calculating_point_of_given_multiplicity([][]int64{Point[0]}, degree)[0]})
		if In(jPs, Q_[0]) {
			l = -(Index(jPs, Q_[0]))
			M := field.p + 1 + 2*m*k + l
			if M == 0 {
				continue
			} else {
				break
			}
		} else {
			Help_point := [][]int64{{Q_[0][0], -(Q_[0][1])}}
			if In(jPs, Help_point[0]) {
				l = Index(jPs, Help_point[0])
				M := field.p + 1 + 2*m*k + l
				if M == 0 {
					continue
				} else {
					break
				}
			}
		}
	}
	M := field.p + 1 + 2*m*k + l
	fact_M := factorize(int64(M), field.p)
	for i := 0; i < len(fact_M); i++ {
		new_M := M / int64(fact_M[i])
		if reflect.DeepEqual(field.Calculating_point_of_given_multiplicity(Point, new_M)[0], []int64{0}) {
			M = new_M
		}
	}
	return M
}

// order returns the order of the elliptic curve, which is the number of points on the curve.
// If the order is not found, it returns the number of points on the curve that were found.
// This function is not guaranteed to terminate and may run indefinitely.
// The function works by randomly generating points on the curve, calculating the LCM of the orders of the points, and then checking if the order divides the LCM.
// If the order divides the LCM, the function returns the order. Otherwise, the function continues to generate points until the order is found or the maximum number of points is reached.
func (field *elliptical_curve) order() int64 {
	var help_i big.Int
	var N int64
	var LCM big.Int
	LCM = *big.NewInt(1)
	flag := false
	var point_counter int64 = 1
	minimum := int64(-((field.p - 1) / 2))
	maximum := int64((field.p - 1) / 2)
	for true {
		x := int64(rand.Intn(int(maximum-minimum+1)) + int(minimum))
		point := field.Point_creator(x)
		if reflect.DeepEqual(point, [][]int64{{0}}) == false {
			if reflect.DeepEqual(point[0], point[1]) {
				point_counter += 1
			} else {
				point_counter += 2
			}
			LCM = *lcm(&LCM, big.NewInt(field.Big_step_for_giant_small_step_for_baby(point)))
			counter := 0
			sqrtP := new(big.Int).Sqrt(big.NewInt(field.p)).Int64()
			for i := field.p + 1 - 2*sqrtP + 1; i < field.p+1+2*sqrtP+1; i++ {
				help_i = *big.NewInt(i)
				if help_i.Mod(&help_i, &LCM).Int64() == 0 {
					N = i
					counter += 1
				}
			}
			if counter == 1 {
				flag = true
				return N
			}
		}
	}
	if flag == false {
		return point_counter
	}
	return 0
}

// Find_simple_group returns a map of all simple subgroups of the elliptic curve. The keys of the map are the orders of the subgroups, and the values are the points of the subgroups.
// The function works by iterating over all possible x-coordinates of the elliptic curve, and for each point, calculating its order and adding it to the map if the order is not already in the map.
// The function then iterates over all the points in the map, and for each point, calculates all the points of the subgroup of that order, and adds them to the map.
// The function returns the map of all simple subgroups of the elliptic curve.
func (field *elliptical_curve) Find_simple_group() map[int64][][]int64 {
	N := field.order()
	Point_0 := [][]int64{{0}}
	list_fact_N := factorize(N, field.p)
	var Point_in_subgroups [][][]int64
	for i := 0; i < len(list_fact_N); i++ {
		Point_in_subgroups = append(Point_in_subgroups, Point_0)
	}
	var subgropus_dict map[int64][][]int64 = make(map[int64][][]int64)
	for i, key := range list_fact_N {
		subgropus_dict[int64(key)] = Point_in_subgroups[i]
	}
	minimum := int64(-((field.p - 1) / 2))
	maximum := int64((field.p - 1) / 2)
	for x := minimum; x < maximum+1; x++ {
		point := field.Point_creator(x)
		if reflect.DeepEqual(point, [][]int64{{0}}) != true {
			point_order := field.Big_step_for_giant_small_step_for_baby(point)
			if In_2(list_fact_N, point_order) && (In(subgropus_dict[point_order], point[0]) == false) && ((In(subgropus_dict[point_order], point[1])) == false) {
				subgropus_dict[point_order] = append(subgropus_dict[point_order], point[0])
				subgropus_dict[point_order] = append(subgropus_dict[point_order], point[1])
				fmt.Println("\nПодгруппа порядка:", point_order)
				for i := int64(0); i < point_order; i++ {
					point_degr := field.Calculating_point_of_given_multiplicity([][]int64{point[0]}, i)
					subgropus_dict[point_order] = append(subgropus_dict[point_order], point_degr[0])
					fmt.Println(point_degr)
				}
			}
		}
	}
	return subgropus_dict
}

// readChoice asks the user to enter a number (1, 2 or 3), reads the input and returns the number.
// If the input is not a number, or not in the range 1-3, it returns an error.
func readChoice() (int, error) {
	reader := bufio.NewReader(os.Stdin)

	fmt.Print("Введите 1, 2 или 3: ")
	input, err := reader.ReadString('\n')
	if err != nil {
		return 0, fmt.Errorf("ошибка чтения ввода: %w", err) // Оборачиваем ошибку для сохранения контекста
	}

	input = strings.TrimSpace(input)

	number, err := strconv.Atoi(input)
	if err != nil {
		return 0, fmt.Errorf("некорректный формат числа: %w", err) // Оборачиваем ошибку
	}

	if number < 1 || number > 3 {
		return 0, errors.New("введено некорректное число. Ожидается 1, 2 или 3") // Возвращаем свою ошибку
	}

	return number, nil
}

// readChoice2 asks the user to enter a number (1 or 2), reads the input and returns the number.
// If the input is not a number, or not in the range 1-2, it returns an error.
func readChoice2() (int, error) {
	reader := bufio.NewReader(os.Stdin)

	fmt.Print("Введите 1 или 2: ")
	input, err := reader.ReadString('\n')
	if err != nil {
		return 0, fmt.Errorf("ошибка чтения ввода: %w", err)
	}

	input = strings.TrimSpace(input)

	number, err := strconv.Atoi(input)
	if err != nil {
		return 0, fmt.Errorf("некорректный формат числа: %w", err)
	}

	if number != 1 && number != 2 {
		return 0, errors.New("введено некорректное число. Ожидается 1 или 2")
	}

	return number, nil
}

// millerRabinTestInt64 tests whether a given int64 number is prime or composite.
// It performs k iterations of the Miller-Rabin primality test.
// If the number is found to be composite, it returns false.
// If the number is probably prime, it returns true.
// The probability of returning a false positive is at most 4^(-k).
func millerRabinTestInt64(n int64, k int) bool {

	if n < 2 {
		return false
	}
	if n == 2 || n == 3 {
		return true
	}
	if n%2 == 0 {
		return false
	}

	nm1 := n - 1
	s := 0
	r := nm1

	for r%2 == 0 {
		r /= 2
		s++
	}

	randSource := rand.NewSource(time.Now().UnixNano())
	rng := rand.New(randSource)

	for i := 0; i < k; i++ {

		a := rng.Int63n(n-3) + 2
		x := powerInt64(a, r, n)

		if x == 1 || x == nm1 {
			continue
		}

		composite := true
		for j := 0; j < s-1; j++ {
			x = (x * x) % n
			if x == 1 {
				return false
			}
			if x == nm1 {
				composite = false
				break
			}
		}

		if composite {
			return false
		}
	}

	return true
}

// powerInt64 calculates the value of (base^exp) % modulus.
// The calculation is done using the exponentiation by squaring algorithm.
// The result is returned as an int64.
func powerInt64(base, exp, modulus int64) int64 {
	result := int64(1)
	base %= modulus
	for exp > 0 {
		if exp%2 == 1 {
			result = (result * base) % modulus
		}
		base = (base * base) % modulus
		exp /= 2
	}
	return result
}

// main provides a command-line interface for interacting with the elliptic curve
// library. The user is prompted to input the values of p, a, and b, and then
// chooses one of three options: constructing the group of points of the
// elliptic curve, calculating a point of given multiplicity, or finding the
// subgroup of a given prime order. The program then prints the result of the
// chosen operation.
func main() {
	var p, a, b int64 // Здесь интерфейс начинается
	fmt.Println("Для любых операций введите значения p, a и b через пробел: ")
	_, err := fmt.Fscan(os.Stdin, &p, &a, &b) // Читаем из stdin и сохраняем в переменные
	if err != nil {
		fmt.Println("Ошибка чтения или преобразования:", err)
		return
	}
	if millerRabinTestInt64(p, 100) == false {
		fmt.Println("p не является простым")
		os.Exit(1)
	}
	Elliptical_curve_work := New_elliptical_curve()
	Elliptical_curve_work.p = p
	Elliptical_curve_work.a = a
	Elliptical_curve_work.b = b
	fmt.Println("Выберите, что вы хотите: \n- Построить группу точек эллиптической кривой или найти её порядок (1); \n- Вычисление точки заданной кратности (2); \n- Нахождение подгрупп простого порядка (3).")
	reader := bufio.NewReader(os.Stdin)
	reader.ReadString('\n')
	choice_1, err_1 := readChoice()
	if err_1 != nil {
		fmt.Println("Ошибка:", err_1)
		return
	}
	if choice_1 == 1 {
		fmt.Println("Вы хотите построить группу (1) или сразу найти её порядок (2)?")
		choice_2, err_2 := readChoice2()
		if err_2 != nil {
			fmt.Println("Ошибка:", err_2)
			return
		}
		if choice_2 == 1 { // Строим группу
			Field_work := New_Field()
			Field_work.p = p
			Field_work.Symmetrical_field_construction()
			ans := Elliptical_curve_work.elliptical_curve_construction(*Field_work)
			if ans > 50 {
				fmt.Printf("Элементы группы: ")
				for i := 0; i <= 50; i++ {
					fmt.Println(Elliptical_curve_work.field_elements[i])
				}
				fmt.Println("Порядок группы:", ans)
			} else {
				fmt.Println("Элементы группы:", Elliptical_curve_work.field_elements)
				fmt.Println("Порядок группы:", ans)
			}
		} else if choice_2 == 2 { // Сразу находим порядок
			fmt.Println("Порядок группы:", Elliptical_curve_work.order())
		}
	} else if choice_1 == 2 {
		var x, y, multiplicity int64
		fmt.Println("Введите x, y и кратность, разделенные пробелами:")
		_, err := fmt.Fscan(os.Stdin, &x, &y, &multiplicity)
		if err != nil {
			fmt.Println("Ошибка чтения или преобразования:", err)
			return
		}
		fmt.Println(Elliptical_curve_work.Calculating_point_of_given_multiplicity([][]int64{{x, y}}, multiplicity))
	} else if choice_1 == 3 {
		Elliptical_curve_work.Find_simple_group()
	}

}
