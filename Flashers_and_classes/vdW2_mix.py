from thermo import *
from thermo.eos import *
from thermo.eos_mix import R2
from fluids.constants import R


class SRKMIX_vdW2_mixing_rule(SRKMIX):

    eos_pure = SRK
    nonstate_constants_specific = ('ms',)
    kwargs_keys = ('kijs', 'lijs',)
    model_id = 10100

    ddelta_dzs = RKMIX.ddelta_dzs
    ddelta_dns = RKMIX.ddelta_dns
    d2delta_dzizjs = RKMIX.d2delta_dzizjs
    d2delta_dninjs = RKMIX.d2delta_dninjs
    d3delta_dninjnks = RKMIX.d3delta_dninjnks

    def __init__(self, Tcs, Pcs, omegas, zs, kijs=None, lijs=None, T=None, P=None, V=None,
                 fugacities=True, only_l=False, only_g=False):
        super().__init__(Tcs=Tcs, Pcs=Pcs, omegas=omegas, zs=zs, kijs=kijs, T=T, P=P, V=V,
                         fugacities=fugacities, only_l=only_l, only_g=only_g)
        self.N = N = len(Tcs)
        cmps = range(N)
        self.Tcs = Tcs
        self.Pcs = Pcs
        self.omegas = omegas
        self.zs = zs
        self.scalar = scalar = type(zs) is list
        b = 0.0
        if kijs is None:
            if scalar:
                kijs = [[0.0] * N for i in cmps]
            else:
                kijs = zeros((N, N))
        self.kijs = kijs
        self.kwargs = {'kijs': kijs}
        if lijs is None:
            if scalar:
                lijs = [[0.0] * N for i in cmps]
            else:
                lijs = zeros((N, N))
        self.lijs = lijs
        self.kwargs = {'lijs': lijs}
        if self.scalar:
            self.ais = [self.c1 * R2 * Tc * Tc / Pc for Tc, Pc in zip(Tcs, Pcs)]
            self.bs = [self.c2 * R * Tc / Pc for Tc, Pc in zip(Tcs, Pcs)]
            ms = [omega * (1.574 - 0.176 * omega) + 0.480 for omega in omegas]
            for i in range(0, len(zs)):
                for j in range(0, len(zs)):
                    b+=self.zs[i] * self.zs[j] * (self.bs[i] + self.bs[j]) * 0.5 * (1 - self.lijs[i][j])
        else:
            Tc_Pc_ratio = Tcs / Pcs
            self.ais = self.c1R2 * Tcs * Tc_Pc_ratio
            self.bs = bs = self.c2R * Tc_Pc_ratio
            ms = omegas * (1.574 - 0.176 * omegas) + 0.480
            for i in range(0, len(self.zs)):
                for j in range(0, len(self.zs)):
                    b += self.zs[i] * self.zs[j] * (self.bs[i] + self.bs[j]) * 0.5 * (1 - self.lijs[i][j])
        self.b = b
        self.ms = ms
        self.delta = self.b

        self.solve(only_l=only_l, only_g=only_g)
        if fugacities:
            self.fugacities()

    def _fast_init_specific(self, other):
        self.ms = other.ms
        b = 0.0
        if self.scalar:
            for i in range(0, len(self.zs)):
                for j in range(0, len(self.zs)):
                    b += self.zs[i] * self.zs[j] * (self.bs[i] + self.bs[j]) * 0.5 * (1 - self.lijs[i][j])
        else:
            for i in range(0, len(self.zs)):
                for j in range(0, len(self.zs)):
                    b += self.zs[i] * self.zs[j] * (self.bs[i] + self.bs[j]) * 0.5 * (1 - self.lijs[i][j])
        self.delta = b

    def bij(self, i, j):
        lij = self.lijs[i][j]
        return (self.bs[i] + self.bs[j]) / 2 * (1 - lij)

    def calculate_bij_matrix(self):
        N = self.N
        bij_matrix = zeros((N, N))
        for i in range(N):
            for j in range(i, N):
                bij_matrix[i][j] = bij_matrix[j][i] = self.bij(i, j)
        return bij_matrix

    def a_alphas_vectorized(self, T):
        return SRK_a_alphas_vectorized(T, self.Tcs, self.aij_matrix, self.ms,
                                       a_alphas=[0.0] * self.N if self.scalar else zeros(self.N))

    def a_alpha_and_derivatives_vectorized(self, T):

        N = self.N
        if self.scalar:
            a_alphas, da_alpha_dTs, d2a_alpha_dT2s = [0.0] * N, [0.0] * N, [0.0] * N
        else:
            a_alphas, da_alpha_dTs, d2a_alpha_dT2s = zeros(N), zeros(N), zeros(N)
        return SRK_a_alpha_and_derivatives_vectorized(T, self.Tcs, self.aij_matrix, self.ms,
                                                      a_alphas=a_alphas, da_alpha_dTs=da_alpha_dTs,
                                                      d2a_alpha_dT2s=d2a_alpha_dT2s)

    def fugacity_coefficients(self, Z):

        N = self.N
        return SRK_lnphis(self.T, self.P, Z, self.b, self.a_alpha, self.bs, self.a_alpha_j_rows, N,
                          lnphis=[0.0] * N if self.scalar else zeros(N))

    def dlnphis_dT(self, phase):

        zs = self.zs
        if phase == 'g':
            Z = self.Z_g
            dZ_dT = self.dZ_dT_g
        else:
            Z = self.Z_l
            dZ_dT = self.dZ_dT_l

        da_alpha_dT_j_rows = self._da_alpha_dT_j_rows
        N = self.N
        P, bs, b = self.P, self.bs, self.b

        T_inv = 1.0 / self.T
        A = self.a_alpha * P * R_inv * T_inv * T_inv
        B = b * P * R_inv * T_inv

        x2 = T_inv * T_inv
        x4 = P * b * R_inv
        x6 = x4 * T_inv

        aij_matrix = self.aij_matrix
        bij_matrix = self.bij_matrix

        da_ij_dT_matrix = zeros((N, N))
        db_ij_dT_matrix = zeros((N, N))

        for i in range(N):
            for j in range(i, N):
                aii = self.ais[i]
                aij = aij_matrix[i][j]
                bij = bij_matrix[i][j]
                bii = self.bs[i]
                db_ij_dT_matrix[i][j] = db_ij_dT_matrix[j][i] = (bii + self.bs[j]) / 2 * da_alpha_dT_j_rows[j]
                if i == j:
                    da_ij_dT_matrix[i][j] = 0
                else:
                    da_ij_dT_matrix[i][j] = da_ij_dT_matrix[j][i] = (aii * aij) ** 0.5 * (
                                da_alpha_dT_j_rows[i] + da_alpha_dT_j_rows[j]) / 2

        db_dT = (db_ij_dT_matrix * zs).sum(axis=1)
        da_dT = (da_ij_dT_matrix * zs).sum(axis=1)

        dlnphis_dT = (aij_matrix / A - 2 * bij_matrix / B) * (da_dT - db_dT)
        return dlnphis_dT
