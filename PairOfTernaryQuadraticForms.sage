"""
@author: Robert Harron

Sage (python) script for working with pairs of ternary quadratic
forms and quartic rings.

"""

load('BinaryCubicForm.sage')
load('TernaryQuadraticForm.sage')

def _matrix_to_TQF_coeffs(M, base_ring=None, integer_matrix=False):
    r"""
    A helper function that converts a matrix to a list of coefficients of a
    ternary quadratic form. If ``integer_matrix`` is ``False``, the matrix
    is interpreted as the Gram matrix of the corresponding quadratic form.
    Otherwise, the (i,j)-entry of the matrix is interpreted as the coefficient
    of x_i*x_j in the quadratic form.

    INPUT:

    - ``M`` -- a 3x3 matrix, assumed to be symmetric
    - ``base_ring`` -- (Default: ``None``) If ``None``, it is set to the
    base ring of ``M``. The entries of ``M`` are turned into coefficients
    of a ternary quadratic form and converted to elements of ``base_ring``.
    - ``integer_matrix`` (Default: False) How to translate the entries of
    ``M`` into coefficients of a ternary quadratic form. If ``True``, the
    (i,j)-entry of the matrix is interpreted as the coefficient of x_i*x_j
    in the quadratic form. If ``False``, the matrix is interpreted as the
    Gram matrix of the ternary quadratic form, so that off-diagonal entries
    are multiplied by 2 when turned into coefficients of the quadratic form.

    OUTPUT:

    A list of length 6 containing the coefficients of a ternary quadratic
    form ordered as [a00, a11, a11, a01, a12, a02].
    """
    if base_ring is None:
        base_ring = M.base_ring()
    coeffs = []
    for i in range(3):
        for j in range(3-i):
            if integer_matrix:
                coeffs.append(base_ring(M[j, j+i]))
            else:
                coeffs.append(base_ring(M[j, j+i] * (2 if i != 0 else 1)))
    return coeffs

def _TQF_coeffs_to_matrix(coeffs, base_ring=None):
    r"""
    A helper function that converts a list of coefficients of a
    ternary quadratic form to the corresonding Gram matrix. The coefficients
    are ordered as [a00, a11, a11, a01, a12, a02].

    INPUT:

    - ``coeffs`` -- list of length 6 ordered as [a00, a11, a11, a01, a12, a02]
    - ``base_ring`` -- (Default: ``None``) If ``None``, it is set to using the
    base ring of the vector [a00, a11, a11, a01/1, a12/1, a02/1]. The base ring
    of the output matrix is set to ``base_ring``.

    OUTPUT:

    A 3x3 symmetrix matrix that is the Gram matrix of the quadratic form.

    ..TODO::

    - Implement ``integer_matrix`` parameter.
    """
    mod_coeffs = [coeffs[i] if i < 3 else coeffs[i] / 2 for i in range(6)]
    if base_ring is None:
        base_ring = vector(mod_coeffs).base_ring()
    A = Matrix(base_ring, 3)
    for i in range(3):
        A[i,i] = mod_coeffs[i]
    for i in range(2):
        A[i, i+1] = mod_coeffs[3+i]
        A[i+1, i] = A[i, i+1]
    A[0,2] = mod_coeffs[-1]
    A[2,0] = A[0,2]
    return A

def _naive_pair_from_poly(f, var_names='x1, x2, x3'):
    r"""
    Helper function that takes in a degree 4 polynomial any outputs
    a pair of ternary quadratic forms that corresponds to the monogenic
    ring the polynomial generates. See the first example for the
    explicit map from polynomial to pair of ternary quadratic forms.

    INPUT:

    - ``f`` -- a degree 4 polynomial (in one variable)
    - ``var_names`` -- (Default: 'x1, x2, x3') the names for the
    variables of the quadratic forms.

    EXAMPLES::

        sage: R.<a0, a1, a2, a3, a4> = PolynomialRing(QQ)
        sage: Rx.<x> = PolynomialRing(R)
        sage: f = Rx(list(R.gens())); f
        a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
        sage: PTQF = _naive_pair_from_poly(f); PTQF.print_pair()
        [     0 1/2*a3 1/2*a1|    -1      0      0]
        [1/2*a3     a4 1/2*a2|     0      0    1/2]
        [1/2*a1 1/2*a2     a0|     0    1/2      0]
    """
    cs = f.list()
    #coeffs = [-1,0,0,0,1,0, cs[2], cs[4], cs[0], cs[3], 0, cs[1]] #this is with y=x^2 and a_2x^2=a_2x^2
    #coeffs = [-1,0,0,0,1,0, 0, cs[0], cs[4], cs[1], cs[2], cs[3]]
    coeffs = [0, cs[4], cs[0], cs[3], cs[2], cs[1], -1,0,0,0,1,0]
    R = vector(coeffs).base_ring()
    return PairOfTernaryQuadraticForms(coeffs, var_names=var_names, base_ring=R)

def _naive_pair_from_poly_flipped(f, var_names='x1, x2, x3'):
    r"""
    A helper function with the same functionality as ``_naive_pair_from_poly``,
    except it switches the order of the output pair.
    """
    cs = f.list()
    coeffs = [-1,0,0,0,1,0, 0, cs[4], cs[0], cs[3], cs[2], cs[1]]
    R = vector(coeffs).base_ring()
    return PairOfTernaryQuadraticForms(coeffs, var_names=var_names, base_ring=R)

# Helper dictionary that converts from (i,j) subscripts to ordering of the coefficients
ij_to_cs = {(1, 1) : 0, (2, 2) : 1, (3, 3) : 2, (1, 2) : 3, (2, 1) : 3, (2, 3) : 4, (3, 2) : 4, (1, 3) : 5, (3, 1) : 5}

class PairOfTernaryQuadraticForms(SageObject):
    r"""
    A class for dealing with pairs of ternary quadratic forms, especially
    as they relate to parameterizing quartic rings. For more on the
    mathematics behind this relation, see Manjul Bhargava's
    'Higher composition laws III: The parametrization of quartic rings'.

    EXAMPLE::

        One way to create a pair of ternary quadratic forms is to pass it
        a list of coefficients. This list has the coefficients of one
        ternary quadratic form and then the other.

        sage: R = PolynomialRing(QQ, names=['a%s%s'%(i, j) for i in range(3) for j in range(3)] + ['b%s%s'%(i, j) for i in range(3) for j in range(3)])
        sage: R.inject_variables()
        Defining a00, a01, a02, a10, a11, a12, a20, a21, a22, b00, b01, b02, b10, b11, b12, b20, b21, b22
        sage: coeffs = [a00, a11, a22, a01, a12, a02, b00, b11, b22, b01, b12, b02]; coeffs
        [a00, a11, a22, a01, a12, a02, b00, b11, b22, b01, b12, b02]
        sage: PairOfTernaryQuadraticForms(coeffs, base_ring=R).print_pair()
        [    a00 1/2*a01 1/2*a02|    b00 1/2*b01 1/2*b02]
        [1/2*a01     a11 1/2*a12|1/2*b01     b11 1/2*b12]
        [1/2*a02 1/2*a12     a22|1/2*b02 1/2*b12     b22]

        The simplest way to create a pair of ternary quadratic forms is
        probably from a quartic polynomial:

        sage: x = polygen(ZZ)
        sage: f = x^4 - x^3 - 3*x^2 + 2*x + 4
        sage: PTQF = _naive_pair_from_poly(f); PTQF
        [-x1*x2 + x2^2 + 2*x1*x3 - 3*x2*x3 + 4*x3^2, -x1^2 + x2*x3]

        This gives a pair of ternary quadratic forms corresponding to the
        quartic ring generated by f. We can see this by working with a
        general quartic. Specifically, we can create a quartic algebra with
        basis the normalized basis of the pair of ternary quadratic forms
        obtained from f, and we can see that the negative of the third basis
        element has characteristic polynomial f:

        sage: R.<a0, a1, a2, a3> = PolynomialRing(QQ)
        sage: Rx.<x> = PolynomialRing(R)
        sage: f2 = Rx(list(R.gens()) + [1]); f2
        x^4 + a3*x^3 + a2*x^2 + a1*x + a0
        sage: PTQF2 = _naive_pair_from_poly(f2); PTQF2
        [a3*x1*x2 + x2^2 + a1*x1*x3 + a2*x2*x3 + a0*x3^2, -x1^2 + x2*x3]
        sage: A = PTQF2.algebra()
        sage: (-A.basis()[2]).characteristic_polynomial()
        x^4 + a3*x^3 + a2*x^2 + a1*x + a0

        We can compare discriminants:

        sage: f.discriminant()
        a1^2*a2^2*a3^2 - 4*a0*a2^3*a3^2 - 4*a1^3*a3^3 + 18*a0*a1*a2*a3^3 - 27*a0^2*a3^4 - 4*a1^2*a2^3 + 16*a0*a2^4 + 18*a1^3*a2*a3 - 80*a0*a1*a2^2*a3 - 6*a0*a1^2*a3^2 + 144*a0^2*a2*a3^2 - 27*a1^4 + 144*a0*a1^2*a2 - 128*a0^2*a2^2 - 192*a0^2*a1*a3 + 256*a0^3
        sage: PTQF2.disc()
        a1^2*a2^2*a3^2 - 4*a0*a2^3*a3^2 - 4*a1^3*a3^3 + 18*a0*a1*a2*a3^3 - 27*a0^2*a3^4 - 4*a1^2*a2^3 + 16*a0*a2^4 + 18*a1^3*a2*a3 - 80*a0*a1*a2^2*a3 - 6*a0*a1^2*a3^2 + 144*a0^2*a2*a3^2 - 27*a1^4 + 144*a0*a1^2*a2 - 128*a0^2*a2^2 - 192*a0^2*a1*a3 + 256*a0^3

        We can also display a pair of ternary quadratic forms as the pair of corresponding Gram matrices, or we can get the matrices separately:

        sage: PTQF2.print_pair()
        [     0 1/2*a3 1/2*a1|    -1      0      0]
        [1/2*a3      1 1/2*a2|     0      0    1/2]
        [1/2*a1 1/2*a2     a0|     0    1/2      0]
        sage: A, B = PTQF2.matrix_pair()
        sage: A
        [     0 1/2*a3 1/2*a1]
        [1/2*a3      1 1/2*a2]
        [1/2*a1 1/2*a2     a0]
        sage: B
        [ -1   0   0]
        [  0   0 1/2]
        [  0 1/2   0]
        (Using PTQF2.show() prints the pair out all pretty!)

        We can ask for a multiplication table:

        sage: PTQF2.multiplication_table()
                        | 1             alpha1   alpha2   alpha3
        +---------------+-------------+--------+--------+--------+
          alpha1^2      | -a1*a3 - a0   -a2      a1       -a3
          alpha1*alpha2 | -a1           0        0        -1
          alpha1*alpha3 | -a0*a3        0        a0       -a2
          alpha2^2      | -a2           -1       a3       0
          alpha2*alpha3 | -a0           0        0        0
          alpha3^2      | 0             a0       0        -a1

        There are several invariants/covariants one has access to
        (for definitions see [HCLIII], or page 114 of Bhargava's
        PhD thesis for the quadratic covariant, which is a
        GL(3)-equivariant ternary quadratic form):

        sage: PTQF2.lambda_ijkl([1,1,2,3])
        a2
        sage: PTQF2.cijk([1,1,3])
        -a3
        sage: PTQF2.cubic_invariant()
        (a1*a2*a3 - a0*a3^2 - a1^2)*X^3 + (-a2^2 - a1*a3 + 4*a0)*X^2*Y + 2*a2*X*Y^2 - Y^3
        sage: PTQF2.quadratic_covariant()
        (4*a2^2 - 8*a1*a3 - 16*a0)*x1^2 + (4*a2*a3 - 24*a1)*x1*x2 + (3*a3^2 - 8*a2)*x2^2 + (4*a1*a2 - 24*a0*a3)*x1*x3 + (2*a1*a3 - 32*a0)*x2*x3 + (3*a1^2 - 8*a0*a2)*x3^2

        We can act by GL(2), GL(3), or a pair (g3, g2) where g3 is in GL(3) and g2 is in GL(2):

        sage: g2 = Matrix(2, [1, 5, 0, 2])
        sage: g3 = Matrix(3, [2, 0, 1, 0, 3, 0, 1, 0, 0])
        sage: PTQF2.act_by_GL2(g2).print_pair()
        [    -6     a3     a1|    -2      0      0]
        [    a3      2 a2 + 3|     0      0      1]
        [    a1 a2 + 3   2*a0|     0      1      0]
        sage: PTQF2.act_by_GL2(g3).print_pair()
        [   0   a3   a1|  -3    0    0]
        [  a3    2   a2|   0    0  3/2]
        [  a1   a2 2*a0|   0  3/2    0]
        sage: PTQF2.action([g3, g2]).print_pair()
        [2*a0 + 4*a1 - 24  3*a2 + 6*a3 + 9          a1 - 12|              -8                3               -4]
        [ 3*a2 + 6*a3 + 9               18             3*a3|               3                0                0]
        [         a1 - 12             3*a3               -6|              -4                0               -2]

        And here's an interesting example of a ring of integers in an S_4 quartic field
        whose cubic resolvent ring is not monogenic:

        sage: PTQF = PairOfTernaryQuadraticForms([0, 2, -2, -2, 1, 2, -1, 1, -1, -1, 1, 1])
        sage: PTQF.disc().factor()
        -1 * 47 * 73
        sage: BCF = PTQF.cubic_invariant(); BCF
        -4*X^3 - 9*X^2*Y + 13*X*Y^2 - 4*Y^3
        sage: [BCF(i, j) for i in GF(2) for j in GF(2)]
        [0, 0, 0, 0]

        Since the discriminant of PTQF is squarefree, we know it corresponds to a maximal
        order. And since BCF doesn't represent 1 mod 2, we know it corresponds to a
        non-monogenic ring!
    """
    def __init__(self, coeffs, base_ring=None, var_names='x1, x2, x3', integer_matrix=False):
        """
        The universal pair:

        EXAMPLES::

            sage: ss = ['a%s%s,'%(i,j) for i in range(1,4) for j in range(i,4)]
            sage: ss += ['b%s%s,'%(i,j) for i in range(1,4) for j in range(i,4)]
            sage: s = ss[0]
            sage: for i in ss[1:]:
                    s += i
            sage: s = s[:-1]
            sage: R = PolynomialRing(QQ, s)
            sage: R.inject_variables()
            sage: AB = PairOfTernaryQuadraticForms([a11,a22,a33,a12,a23,a13,b11,b22,b33,b12,b23,b13],var_names='x,y,z')
            sage: AB
            [a11*x^2 + a12*x*y + a22*y^2 + a13*x*z + a23*y*z + a33*z^2, b11*x^2 + b12*x*y + b22*y^2 + b13*x*z + b23*y*z + b33*z^2]
        """
        if len(coeffs) == 2:
            self._as = _matrix_to_TQF_coeffs(coeffs[0], integer_matrix=integer_matrix)
            self._bs = _matrix_to_TQF_coeffs(coeffs[1], integer_matrix=integer_matrix)
            if base_ring is None:
                base_ring = vector(self._as + self._bs).base_ring()
        else:
            if base_ring is None:
                base_ring = vector(coeffs).base_ring()
            self._as = coeffs[:6]    #stored as a11, a22, a33, a12, a23, a13, the coefficients of xixj
            self._bs = coeffs[6:12]
        self._as = [base_ring(a) for a in self._as]
        self._bs = [base_ring(b) for b in self._bs]
        self._base_ring = base_ring
        self._polynomial_ring = PolynomialRing(base_ring, var_names)
        self._vars = self._polynomial_ring.gens()
        self._polynomials = [sum([self._as[(-i^2 / 2 + 7 * i / 2) + j] * self._vars[j] * self._vars[j+i] for i in range(3) for j in range(3 - i)]), sum([self._bs[(-i^2 / 2 + 7 * i / 2) + j] * self._vars[j] * self._vars[j+i] for i in range(3) for j in range(3 - i)])] #is this different for integer_matrix?
        self._int_mat = integer_matrix
        return

    def _repr_(self):
        r"""
        Returns the string representation of list of polynomials underlying ``self``.
        """
        return str(self._polynomials)

    def show(self):
        r"""
        Displays a 'pretty' output of the pair of Gram matrices corresponding to ``self``.
        """
        A, B = self.matrix_pair()
        show(A,B)
        return

    def coefficients(self):
        r"""
        Returns a list of length 12 containg the coefficients of ``self``.

        EXAMPLES::

            sage: PairOfTernaryQuadraticForms([1,2,3,4,5,6,7,8,9,10,11,12]).coefficients()
            [1,2,3,4,5,6,7,8,9,10,11,12]
        """
        return self._as + self._bs

    def matrix_pair(self):
        r"""
        Returns the pair of Gram matrices (as a tuple) corresponding to the two
        ternary quadratic forms underlying ``self``.

        This function attempts to preserve the base ring of ``self``, but might
        not be able to because of denominators of 2 in the Gram matrix.
        """
        try:
            A = _TQF_coeffs_to_matrix(self._as, base_ring=self._base_ring)
            B = _TQF_coeffs_to_matrix(self._bs, base_ring=self._base_ring)
        except TypeError:
            A = _TQF_coeffs_to_matrix(self._as)
            B = _TQF_coeffs_to_matrix(self._bs)
        return (A, B)

    def print_pair(self):
        r"""
        Prints out the pair of Gram matrices corresponding to the two
        ternary quadratic forms underlying ``self``. They are displayed
        next to each other.

        EXAMPLES::

            sage: PairOfTernaryQuadraticForms([1,2,3,4,5,6,7,8,9,10,11,12]).print_pair()
            [   1    2    3|   7    5    6]
            [   2    2  5/2|   5    8 11/2]
            [   3  5/2    3|   6 11/2    9]
        """
        A, B = self.matrix_pair()
        print A.augment(B, subdivide=True)
        return

    def lambda_ijkl(self, ijkl):
        r"""
        Return the lambda^{ij}_{kl} invariants of ``self`` as defined in equation (20)
        of [HCL III].

        INPUT:

        - ``ijkl`` -- a list of length 4 containing integers between 1 and 3.
        """
        i, j, k, l = ijkl
        ij = (i, j)
        kl = (k, l)
        return Matrix(2, [self._as[ij_to_cs[ij]], self._bs[ij_to_cs[ij]], self._as[ij_to_cs[kl]], self._bs[ij_to_cs[kl]]]).det()

    def lambdas(self):
        r"""
        Return the list of all lambda^{ij}_{kl} invariants as a dictionary whose
        keys are the tuples (i,j,k,l).

        Note that there are only 15 lambda-invariants since certain permutations
        of {i,j,k,l} give the same invariants. The output of this function is a
        dictionary with 15 entries whose keys are
        (1, 1, 1, 2)
        (1, 1, 1, 3)
        (1, 1, 2, 2)
        (1, 1, 2, 3)
        (1, 1, 3, 3)
        (1, 2, 1, 3)
        (1, 2, 2, 2)
        (1, 2, 2, 3)
        (1, 2, 3, 3)
        (1, 3, 2, 2)
        (1, 3, 2, 3)
        (1, 3, 3, 3)
        (2, 2, 2, 3)
        (2, 2, 3, 3)
        (2, 3, 3, 3)
        """
        Ls = []
        for i in range(1,4):
            for j in range(i,4):
                for k in range(i, 4):
                    for l in range(k,4):
                        if k == i and j >= l:
                            continue
                        ijkl = (i, j, k, l)
                        Ls.append([ijkl, self.lambda_ijkl(ijkl)])
        return dict(Ls)

    def cijk(self, ijk):
        r"""
        Return the c^k_{ij} invariants of ``self`` defined in equation (14)
        of [HCL III] as the coefficients arising in the multiplication table
        of the normalized basis of the quartic ring corresponding to ``self``.

        INPUT:

        - ``ijk`` -- a triple of integers between 0 and 3.
        """
        i, j, k = ijk
        ijk = [1,2,3]
        if k == 0:
            if i > j:
                temp = i
                i = j
                j = temp
            ijk.remove(j)
            if i != j:
                ijk.remove(i)
            k = ijk[0]
            return sum([self.cijk([j,k,r]) * self.cijk([r,i,k]) - self.cijk([i,j,r]) * self.cijk([r,k,k]) for r in range(1, 4)])
        if j == k and i != j:
            j = i
            i = k
        if i == j:
            if j == k: #ciii
                ijk.remove(i)
                j, k = ijk
                return Permutation([i, j, k]).sign() * self.lambda_ijkl([i,k,i,j]) + self.Ci(i)
            else: #ciij
                j = k
                ijk.remove(i)
                ijk.remove(j)
                k = ijk[0]
                return Permutation([i, j, k]).sign() * self.lambda_ijkl([i,i,i,k])
        else:
            if i == k: #ciji
                ijk.remove(i)
                ijk.remove(j)
                k = ijk[0]
                return Permutation([i, j, k]).sign() * self.lambda_ijkl([i,k,j,j]) / 2 + self.Ci(j) / 2
            else: #cijk
                return Permutation([i, j, k]).sign() * self.lambda_ijkl([j,j,i,i])

    def Ci(self, i):
        r"""
        Return the C_i invariant given in equation (22) of [HCL III].

        INPUT:

        - ``i`` -- an integer between 1 and 3.
        """
        if i == 1:
            return self.lambda_ijkl([2,3,1,1])
        if i == 2:
            return -self.lambda_ijkl([1,3,2,2])
        return self.lambda_ijkl([1,2,3,3])

    def multiplication_table(self, names='alpha'):
        r"""
        Return the multiplication table of the normalized basis of
        the quartic ring corresponding to ``self``.

        INPUT:

        - ``names`` -- (Default: 'alpha') the basis elements will be
        labelled using ``names`` with subscripts 0 to 3.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: f = x^4 - 2
            sage: PTQF = _naive_pair_from_poly(f)
            sage: PTQF.multiplication_table()
                            | 1   alpha1   alpha2   alpha3
            +---------------+---+--------+--------+--------+
              alpha1^2      | 2   0        0        0
              alpha1*alpha2 | 0   0        0        -1
              alpha1*alpha3 | 0   0        -2       0
              alpha2^2      | 0   -1       0        0
              alpha2*alpha3 | 2   0        0        0
              alpha3^2      | 0   -2       0        0
        """
        names2 = [names + '%s'%(i) for i in range(1,4)]
        R = PolynomialRing(QQ, names2)
        R.inject_variables(verbose=False);
        return table([['', '1'] + names2] + [[eval(names + '%s*'%(i) + names +'%s'%(j))] + [self.cijk([i,j,k]) for k in range(4)] for i in range(1,4) for j in range(i,4)], header_row=True, header_column=True)

    def algebra(self, names='alpha'):
        r"""
        Returns a FiniteDimensionalAlgebra over the fraction field of ``self.base_ring()``
        (so the latter must be an integral domain). This algebra has basis the normalized
        basis of the quartic ring corresponding to ``self``. This gives access to all the
        functionality of FiniteDimensionalAlgebra.

        INPUT:

        - ``names`` -- (Default: 'alpha') the basis elements of the outputted algebra will
        be labelled using ``names`` with subscripts 0 to 3.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: f = x^4 - 2
            sage: PTQF = _naive_pair_from_poly(f)
            sage: A = PTQF.algebra()
            sage: for alpha in A.basis():
                      alpha.minimal_polynomial()
            x - 1
            x^2 - 2
            x^4 - 2
            x^4 - 8
            sage: (A.basis()[1] + A.basis()[2]).matrix().det()
            2

            The latter computation gives the norm of 2^(1/2) + 2^(1/4).
        """
        matrices = [identity_matrix(4)] + [matrix(4, [self.cijk([i,j,k]) if j != 0 else ZZ(i == k) for k in range(4) for j in range(4)]).transpose() for i in range(1,4)]
        return FiniteDimensionalAlgebra(self._base_ring.fraction_field(), matrices, assume_associative=True, names=names)

    def cubic_invariant(self, var_names='X, Y'):
        r"""
        Returns the cubic GL(2)-invariant of ``self``.

        Letting A and B denote the Gram matrices of the ternary
        quadratic forms underlying ``self``, the output is a
        BinaryCubicForm given by

        4 * det(A*X - B*Y)

        (with the factor of 4 removed if ``self`` was created with
        integer_matrix = True).

        INTPUT:

        - ``var_names`` -- (Default: 'X, Y') the variable names to
        use for the outputted binary cubic form.
        """
        A, B = self.matrix_pair()
        _.<X, Y> = PolynomialRing(A.base_ring(), var_names)
        M = A*X - B*Y
        D = M.det()
        if not self._int_mat:
            D *= 4
        coeffs = [D[X^(3-i)*Y^i] for i in range(4)]
        return BinaryCubicForm(coeffs, base_ring=self._base_ring, var_names=var_names)

    def quadratic_covariant(self, var_names=None):
        r"""
        Returns the quadratic GL(3)-covariant of ``self``.

        A formula for this is given on page 114 of Bhargava's PhD thesis. It
        is 4 times the Gram matrix of the trace-zero form of the quartic
        ring corresponding to ``self`` with respect to its normalized basis.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: f = x^4 - 2
            sage: PTQF = _naive_pair_from_poly(f)
            sage: K.<a2> = NumberField(f)
            sage: a1 = -a2^2
            sage: a3 = -a1 * a2
            sage: K.trace_pairing([1, a1, a2, a3])
            [4 0 0 0]
            [0 8 0 0]
            [0 0 0 8]
            [0 0 8 0]
            sage: PTQF.quadratic_covariant()
            32*x1^2 + 64*x2*x3
        """
        if var_names is None:
            var_names = ','.join(self._polynomial_ring.variable_names())
        lambda_ijkl =  self.lambda_ijkl
        Q11 = (-16 * lambda_ijkl([1,1,2,2])*lambda_ijkl([1,1,3,3]) - 4 * lambda_ijkl([1,1,1,3])*lambda_ijkl([1,2,2,3])
               + 8 * lambda_ijkl([1,1,1,3])*lambda_ijkl([1,3,2,2]) - 4 * lambda_ijkl([1,1,1,2]) * lambda_ijkl([1,3,2,3])
               + 8 * lambda_ijkl([1,1,1,2])*lambda_ijkl([1,2,3,3]) + 4 * lambda_ijkl([1,1,2,3])^2 + 3 * lambda_ijkl([1,2,1,3])^2)
        Q22 = (16 * lambda_ijkl([1,1,2,2])*lambda_ijkl([2,2,3,3]) - 4 * lambda_ijkl([1,2,1,3])*lambda_ijkl([2,2,2,3])
               - 8 * lambda_ijkl([1,1,2,3])*lambda_ijkl([2,2,2,3]) - 4 * lambda_ijkl([1,2,2,2]) * lambda_ijkl([1,3,2,3])
               - 8 * lambda_ijkl([1,2,2,2])*lambda_ijkl([1,2,3,3]) + 4 * lambda_ijkl([1,3,2,2])^2 + 3 * lambda_ijkl([1,2,2,3])^2)
        Q33 = (-16 * lambda_ijkl([1,1,3,3])*lambda_ijkl([2,2,3,3]) - 4 * lambda_ijkl([1,2,1,3])*lambda_ijkl([2,3,3,3])
               + 8 * lambda_ijkl([1,1,2,3])*lambda_ijkl([2,3,3,3]) - 4 * lambda_ijkl([1,2,2,3]) * lambda_ijkl([1,3,3,3])
               - 8 * lambda_ijkl([1,3,2,2])*lambda_ijkl([1,3,3,3]) + 4 * lambda_ijkl([1,2,3,3])^2 + 3 * lambda_ijkl([1,3,2,3])^2)
        Q12 = (8 * lambda_ijkl([1,1,1,2])*lambda_ijkl([2,2,3,3]) - 8 * lambda_ijkl([1,1,3,3])*lambda_ijkl([1,2,2,2])
               -12 * lambda_ijkl([1,1,1,3])*lambda_ijkl([2,2,2,3]) + 2 * lambda_ijkl([1,1,2,3]) * lambda_ijkl([1,2,2,3])
               - 2 * lambda_ijkl([1,2,1,3])*lambda_ijkl([1,3,2,2]) + 1 * lambda_ijkl([1,2,1,3]) * lambda_ijkl([1,2,2,3]))
        Q23 = (8 * lambda_ijkl([1,1,2,2])*lambda_ijkl([2,3,3,3]) - 8 * lambda_ijkl([1,1,3,3]) * lambda_ijkl([2,2,2,3])
               -12 * lambda_ijkl([1,2,2,2])*lambda_ijkl([1,3,3,3]) + 2 * lambda_ijkl([1,2,2,3]) * lambda_ijkl([1,2,3,3])
               + 2 * lambda_ijkl([1,3,2,2])*lambda_ijkl([1,3,2,3]) + 1 * lambda_ijkl([1,2,2,3]) * lambda_ijkl([1,3,2,3]))
        Q13 = (-8 * lambda_ijkl([1,1,1,3])*lambda_ijkl([2,2,3,3]) - 8 * lambda_ijkl([1,1,2,2])*lambda_ijkl([1,3,3,3])
               +12 * lambda_ijkl([1,1,1,2])*lambda_ijkl([2,3,3,3]) + 2 * lambda_ijkl([1,1,2,3]) * lambda_ijkl([1,3,2,3])
               + 2 * lambda_ijkl([1,2,1,3])*lambda_ijkl([1,2,3,3]) - 1 * lambda_ijkl([1,2,1,3]) * lambda_ijkl([1,3,2,3]))
        return TernaryQuadraticForm([Q11, Q22, Q33, 2*Q12, 2*Q23, 2*Q13], base_ring=self._base_ring, var_names=','.join(self._polynomial_ring.variable_names()))

    def disc(self):
        r"""
        Return the discriminant of ``self``.

        By definition, this is the discriminants of its cubic invariant. It
        also equals the discriminant of the corresponding quartic ring.
        """
        return self.cubic_invariant().disc()

    def act_by_GL2(self, gamma):
        r"""
        Return the PairOfTernaryQuadraticForms obtained by acting on
        ``self`` by the 2x2 invertible matrix ``gamma``.
        """
        a2 = gamma[0,0] * vector(self._as) + gamma[0,1] * vector(self._bs)
        b2 = gamma[1,0] * vector(self._as) + gamma[1,1] * vector(self._bs)
        return PairOfTernaryQuadraticForms(list(a2) + list(b2), var_names=','.join(self._polynomial_ring.variable_names()))

    def act_by_GL3(self, gamma):
        r"""
        Return the PairOfTernaryQuadraticForms obtained by acting on
        ``self`` by the 3x3 invertible matrix ``gamma``.
        """
        A, B = self.matrix_pair()
        gAgt = gamma * A * gamma.transpose()
        gBgt = gamma * B * gamma.transpose()
        try:
            a2 = _matrix_to_TQF_coeffs(gAgt, base_ring=self._base_ring)
            b2 = _matrix_to_TQF_coeffs(gBgt, base_ring=self._base_ring)
        except TypeError:
            a2 = _matrix_to_TQF_coeffs(gAgt)
            b2 = _matrix_to_TQF_coeffs(gBgt)
        return PairOfTernaryQuadraticForms(a2 + b2, var_names=','.join(self._polynomial_ring.variable_names()))

    def action(self, gs):
        r"""
        Return the PairOfTernaryQuadraticForms obtained by acting on
        ``self`` by the pair (g3, g2) where g3 is a 3x3 invertiable
        matrix and g2 is a2x2 invertible matrix.

        This is the same as first acting by one of them and then the other.
        """
        g3, g2 = gs
        return self.act_by_GL2(g2).act_by_GL3(g3)
