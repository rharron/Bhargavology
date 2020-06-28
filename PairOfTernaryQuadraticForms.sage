"""
@author: Robert Harron

Sage (python) script for working with pairs of ternary quadratic
forms and quartic rings.

"""

load('BinaryCubicForm.sage')
load('TernaryQuadraticForm.sage')

def _matrix_to_TQF_coeffs(M, base_ring=None, integer_matrix=False):
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

def _TQF_coeffs_to_matrix(coeffs, base_ring=None):#, integer_matrix=False):
    #if not integer_matrix:
    #    mod_coeffs = [coeffs[i] if i < 3 else coeffs[i] / 2 for i in range(6)]
    #else:
    #    mod_coeffs = coeffs
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
    cs = f.list()
    #coeffs = [-1,0,0,0,1,0, cs[2], cs[4], cs[0], cs[3], 0, cs[1]] #this is with y=x^2 and a_2x^2=a_2x^2
    #coeffs = [-1,0,0,0,1,0, 0, cs[0], cs[4], cs[1], cs[2], cs[3]]
    coeffs = [0, cs[4], cs[0], cs[3], cs[2], cs[1], -1,0,0,0,1,0]
    R = vector(coeffs).base_ring()
    return PairOfTernaryQuadraticForms(coeffs, var_names=var_names, base_ring=R)

def _naive_pair_from_poly_flipped(f, var_names='x1, x2, x3'):
    cs = f.list()
    #coeffs = [-1,0,0,0,1,0, cs[2], cs[4], cs[0], cs[3], 0, cs[1]] #this is with y=x^2 and a_2x^2=a_2x^2
    #coeffs = [-1,0,0,0,1,0, 0, cs[0], cs[4], cs[1], cs[2], cs[3]]
    coeffs = [-1,0,0,0,1,0, 0, cs[4], cs[0], cs[3], cs[2], cs[1]]
    R = vector(coeffs).base_ring()
    return PairOfTernaryQuadraticForms(coeffs, var_names=var_names, base_ring=R)

ij_to_cs = {(1, 1) : 0, (2, 2) : 1, (3, 3) : 2, (1, 2) : 3, (2, 1) : 3, (2, 3) : 4, (3, 2) : 4, (1, 3) : 5, (3, 1) : 5}

class PairOfTernaryQuadraticForms(SageObject):
    r"""
    A class for dealing with pairs of ternary quadratic forms, especially
    as they relate to parameterizing quartic rings.
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
        return str(self._polynomials)

    def show(self):
        A, B = PTQF.matrix_pair()
        show(A,B)
        return

    def coefficients(self):
        return self._as + self._bs

    def matrix_pair(self):
        try:
            A = _TQF_coeffs_to_matrix(self._as, base_ring=self._base_ring)
            B = _TQF_coeffs_to_matrix(self._bs, base_ring=self._base_ring)
        except TypeError:
            A = _TQF_coeffs_to_matrix(self._as)
            B = _TQF_coeffs_to_matrix(self._bs)
        return (A, B)

    def print_pair(self):
        A, B = self.matrix_pair()
        print A.augment(B, subdivide=True)
        return

    def lambda_ijkl(self, ijkl):
        i, j, k, l = ijkl
        ij = (i, j)
        kl = (k, l)
        return Matrix(2, [self._as[ij_to_cs[ij]], self._bs[ij_to_cs[ij]], self._as[ij_to_cs[kl]], self._bs[ij_to_cs[kl]]]).det()

    def lambdas(self):
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

    def hij(self, ij):
        ij = copy(ij)
        ij.sort()
        lambda_ijkl = self.lambda_ijkl
        if ij == [1, 1]:
            return -lambda_ijkl([1,1,2,2])*lambda_ijkl([1,1,3,3])-lambda_ijkl([1,1,1,3])*lambda_ijkl([1,2,2,3])\
             +lambda_ijkl([1,1,1,3])*lambda_ijkl([1,3,2,2]) + lambda_ijkl([1,1,1,2])*lambda_ijkl([1,2,3,3])
        elif ij == [2,2]:
            return lambda_ijkl([1,1,2,2])*lambda_ijkl([2,2,3,3])-lambda_ijkl([1,2,1,3])*lambda_ijkl([2,2,2,3])\
             -lambda_ijkl([1,1,2,3])*lambda_ijkl([2,2,2,3])
        elif ij == [3,3]:
            return -lambda_ijkl([1,1,3,3])*lambda_ijkl([2,2,3,3])-lambda_ijkl([1,2,1,3])*lambda_ijkl([2,3,3,3])
        elif ij == [1,2]:
            return -lambda_ijkl([1,1,1,3])*lambda_ijkl([2,2,2,3])+lambda_ijkl([1,1,1,2])*lambda_ijkl([2,2,3,3])
        elif ij == [1,3]:
            return -lambda_ijkl([1,1,1,3])*lambda_ijkl([2,2,3,3])+lambda_ijkl([1,1,1,2])*lambda_ijkl([2,3,3,3])
        else:
            return -lambda_ijkl([1,1,3,3])*lambda_ijkl([2,2,2,3])-lambda_ijkl([1,2,1,3])*lambda_ijkl([2,2,3,3])

    def Ci(self, i):
        if i == 1:
            return self.lambda_ijkl([2,3,1,1])
        if i == 2:
            return -self.lambda_ijkl([1,3,2,2])
        return self.lambda_ijkl([1,2,3,3])

    def multiplication_table(self, names='alpha'):
        names2 = [names + '%s'%(i) for i in range(1,4)]
        R = PolynomialRing(QQ, names2)
        R.inject_variables(verbose=False);
        return table([['', '1'] + names2] + [[eval(names + '%s*'%(i) + names +'%s'%(j))] + [self.cijk([i,j,k]) for k in range(4)] for i in range(1,4) for j in range(i,4)], header_row=True, header_column=True)

    def algebra(self, names='alpha'):
        matrices = [identity_matrix(4)] + [matrix(4, [self.cijk([i,j,k]) if j != 0 else ZZ(i == k) for k in range(4) for j in range(4)]).transpose() for i in range(1,4)]
        return FiniteDimensionalAlgebra(self._base_ring.fraction_field(), matrices, assume_associative=True, names=names)

    def cubic_invariant(self, var_names='X, Y'):
        A, B = self.matrix_pair()
        _.<X, Y> = PolynomialRing(A.base_ring(), var_names)
        M = A*X - B*Y
        D = M.det()
        if not self._int_mat:
            D *= 4
        coeffs = [D[X^(3-i)*Y^i] for i in range(4)]
        return BinaryCubicForm(coeffs, base_ring=self._base_ring, var_names=var_names)

    def quadratic_covariant(self, var_names=None):
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
        return self.cubic_invariant().disc()

    def act_by_GL2(self, gamma):
        a2 = gamma[0,0] * vector(self._as) + gamma[0,1] * vector(self._bs)
        b2 = gamma[1,0] * vector(self._as) + gamma[1,1] * vector(self._bs)
        return PairOfTernaryQuadraticForms(list(a2) + list(b2), var_names=','.join(self._polynomial_ring.variable_names()))

    def act_by_GL3(self, gamma):
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
        g3, g2 = gs
        return self.act_by_GL2(g2).act_by_GL3(g3)
