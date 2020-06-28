"""
@author: Robert Harron

Sage (python) script for working with ternary quadratic forms
as they relate to pairs of ternary quadratic forms.

"""

class TernaryQuadraticForm(SageObject):
    r"""
    A class for computing with ternary quadratic forms as they
    relate to the code from PairOfTernaryQuadraticForms.

    ..NOTE::This is currently a very bare bones implementation.
    """
    def __init__(self, coeffs, base_ring=None, var_names='x1, x2, x3'):
        if base_ring is None:
            base_ring = vector(coeffs).base_ring()
        self.coeffs = [base_ring(a) for a in coeffs]
        self._base_ring = base_ring
        self._polynomial_ring = PolynomialRing(base_ring, var_names)
        self._vars = self._polynomial_ring.gens()
        self._polynomial = sum([self.coeffs[(-i^2 / 2 + 7 * i / 2) + j] * self._vars[j] * self._vars[j+i] for i in range(3) for j in range(3 - i)])
        return

    def _repr_(self):
        return str(self._polynomial)

    def coefficients(self):
        return self.coeffs

    def matrix(self):
        try:
            A = _TQF_coeffs_to_matrix(self.coeffs, base_ring=self._base_ring)
        except TypeError:
            A = _TQF_coeffs_to_matrix(self.coeffs)
        return A

    def disc(self):
        return 4 * self.matrix().det()

    def action(self, gamma):
        A = self.matrix()
        gAgt = gamma * A * gamma.transpose()
        #gBgt = gamma * B * gamma.transpose()
        try:
            a2 = _matrix_to_TQF_coeffs(gAgt, base_ring=self._base_ring)
            #b2 = _matrix_to_TQF_coeffs(gBgt, base_ring=self._base_ring)
        except TypeError:
            a2 = _matrix_to_TQF_coeffs(gAgt)
            #b2 = _matrix_to_TQF_coeffs(gBgt)
        return TernaryQuadraticForm(a2, var_names=','.join(self._polynomial_ring.variable_names()))

    def primitive_form(self):
        g = gcd(self.coeffs)
        return TernaryQuadraticForm([a/g for a in self.coeffs], base_ring=self._base_ring, var_names=','.join(self._polynomial_ring.variable_names()))
