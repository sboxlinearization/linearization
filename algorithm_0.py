import random
import itertools
from collections import defaultdict, Counter

from binteger import Bin

from affine import LinearSpan, AffineSpan, XOR

from tqdm import tqdm
from math import factorial

try:
    from sage.all import ZZ, binomial, QQ, copy
    from sage.crypto.sboxes import SBox
except ImportError:
    pass

class Algorithm:
    S: list[int] # the S-box
    n: int # input bits
    m: int # output bits
    def __init__(self, sbox, n, m):
        self.S = tuple(int(y) for y in sbox)
        self.n = int(n)
        self.m = int(m)
        if len(self.S) != 2**self.n:
            raise ValueError("size of the sbox input does not match")
        if not all(0 <= y < 2**self.m for y in self.S):
            raise ValueError("size of the sbox output does not match")

        self.kddt = {}
        self.domain = tuple(range(2**self.n))
        self.codomain = tuple(range(2**self.m))
        self.graph = tuple( (x << self.m) | self.S[x] for x in self.domain )
        self.graph_set = set(self.graph)

        #self.S_inv = None
        #if self.n == self.m and set(self.S) == set(range(2**self.n)):
        #    self.S_inv = [None] * 2**self.m
        #    for x in self.domain:
        #        self.S_inv[self.S[x]] = x

        #self.TAB_SCALAR_X = [[0] * 2**self.n for _ in range(2**self.n)]
        #for a in range(2**self.n):
        #    for b in range(2**self.n):
        #        self.TAB_SCALAR_X[a][b] = bin(a & b).count("1") & 1

        #self.TAB_SCALAR_Y = [[0] * 2**self.m for _ in range(2**self.m)]
        #for a in range(2**self.m):
        #    for b in range(2**self.m):
        #        self.TAB_SCALAR_Y[a][b] = bin(a & b).count("1") & 1

        self._lat = None
        self._ddt = None

    def _make_parity_tab(self, nbits):
        return tuple(bin(x).count("1") & 1 for x in range(2**nbits))

    def lat(self):
        if self._lat is None:
            self._lat = SBox(self.S).linear_approximation_table().change_ring(ZZ)  # bug in Sage: it is QQ by default
        return self._lat

    def ddt(self):
        if self._ddt is None:
            self._ddt = SBox(self.S).difference_distribution_table()
        return self._ddt

    def to_graph(self, xs):
        return tuple( (x << self.m) | self.S[x] for x in xs )

    def is_affine_approximation(self, xs, check_maximal=True):
        asp = AffineSpan(xs)
        aspg = set(AffineSpan(self.to_graph(xs)))
        if len(aspg) != len(asp):
            return False
        if check_maximal:
            if len(asp) != 2**(self.n) or len(aspg & self.graph_set) != len(xs):
                return False
        return True

    def expand_approximation(self, xs):
        asp = set(AffineSpan(self.to_graph(xs)))
        if len(asp) > 2**self.n:
            raise ValueError("inconsistent approximation")
        A_graph = asp & self.graph_set
        return tuple(sorted(xy >> self.m for xy in A_graph))

    def approximation_from_inputs(self, xs):
        asp = AffineSpan(self.to_graph(xs))
        if len(asp) != 2**self.n:
            raise ValueError("approximation is not maximal")
        A = [None] * 2**self.n
        mask = 2**self.m - 1
        for xy in asp:
            x = xy >> self.m
            y = xy & mask
            A[x] = y
        return A

    def sample_kzerosum(self, k=4, affine_basis=True):
        if affine_basis:
            if k % 2 == 0 and 2**(k-2) > 2**self.n:
                raise ValueError("too small S-box for %d-zerosum" % k)
            if k % 2 == 1 and 2**(k-1) > 2**self.n:
                raise ValueError("too small S-box for %d-zerosum" % k)
        if k > 2**self.n:
            raise ValueError("too small S-box for %d-zerosum" % k)
        while True:
            xs = random.sample(self.domain, k-1)
            if affine_basis and len(AffineSpan(xs)) < 2**(k-2):
                continue
            xsum = XOR(xs)
            if affine_basis == False and xsum in set(xs):
                continue
            ysum = XOR(self.S[x] for x in xs)
            if self.S[xsum] == ysum:
                xs.append(xsum)
                xs.sort()
                if k % 2 == 1 and len(AffineSpan(xs)) != 2**(k-1):
                    continue
                return tuple(xs)

    def compute_kddt(self, k=2):
        if k not in self.kddt:
            print("ddt_empty_start")
            tab = [
                #[list() for _ in range(2**self.m)]
                defaultdict(list)
                for _ in tqdm(range(2**self.n))
            ]
            print("ddt_empty_stop")
            # IMPORTANT: lists and tuples should be sorted
            print("ddt_start")
            for xs in tqdm(itertools.combinations(self.domain, k),total=len(self.domain)**k/factorial(k)):
                tab[XOR(xs)][XOR(self.S[x] for x in xs)].append(xs)
            self.kddt[k] = tab
        return self.kddt[k]

    def iter_kzerosums(self, k):
        if k % 2:
            raise NotImplementedError()
        ddt = self.compute_kddt(k=k//2)
        # IMPORTANT: assumes that ddt entries (lists) are sorted,
        # as well as the tuples inside lists
        for row in ddt:
            for lst in row:
                # enumerate pair of tuples (a,...,b), (c,...,d)
                # such that a < ... < b < c < ... < d
                i = 0
                while i < len(lst) - 1:
                    j = i + 1
                    while j < len(lst) and lst[j][0] <= lst[i][-1]:
                        j += 1
                    for jj in range(j, len(lst)):
                        yield lst[i] + lst[jj]
                    i += 1

    def _lat_xors_max_output(self):
        """Old basic version of the algorithm"""
        lat = self.lat()

        res = []
        for cx in self.domain:
            cur_lat = copy(lat)

            # flip rows based on the input constant cx
            tab_scalar_cx = self.TAB_SCALAR_X[cx]
            for mx in self.domain:
                if tab_scalar_cx[mx]:
                    cur_lat[mx] = -cur_lat[mx]

            maxes = []
            for my in self.codomain:
                col = cur_lat.column(my)
                max_pos = max(col)
                max_neg = max(-col)
                maxes.append((max_pos, max_neg))

            for cy in self.codomain:
                cur = 0
                tab_scalar_cy = self.TAB_SCALAR_Y[cy]
                for my in range(1, 2**self.m):
                    cur += maxes[my][tab_scalar_cy[my]]

                res.append((cur, cx, cy))

        res.sort(reverse=True)
        return res

    def lat_xors(self, score_outputs=True, score_inputs=True, exclude_cx_cy=False):
        """Return dict/Counter of {(cx, cy): score}.

        (cx, cy) refer to i/o constants leading to 0 mapping to 0.
        Equivalently, this is an i/o pair candidate of A: A(0+cx) + cy = 0 .
        Score is the sum of best correlations of this XOR over all INPUT+OUTPUT components (lin. comb.),
        depending on the arguments.
        """
        lat = self.lat()
        result = Counter()
        if score_outputs:
            self._lat_xors_score_outputs(lat, result=result, exclude_cx_cy=exclude_cx_cy)
        if score_inputs:
            self._lat_xors_score_outputs(lat.transpose(), result=result, flip_input_output=True, exclude_cx_cy=exclude_cx_cy)
        return result

    def _lat_xors_score_outputs(self, lat, result=None, flip_input_output=False, exclude_cx_cy=False):
        """Return dict/Counter of {(cx, cy): score}.

        (cx, cy) refer to i/o constants leading to 0 mapping to 0.
        Equivalently, this is an i/o pair candidate of A: A(0+cx) + cy = 0 .
        Score is the sum of best correlations of this XOR over all OUTPUT components (lin. comb.).
        """
        n = int(lat.nrows()).bit_length() - 1
        m = int(lat.ncols()).bit_length() - 1
        assert lat.nrows() == 2**n
        assert lat.ncols() == 2**m

        parity_mx = self._make_parity_tab(n)
        parity_my = self._make_parity_tab(m)

        if result is None:
            result = Counter()
        columns = []
        for my in range(2**m):
            # sorted column of (mx, lat[mx,my])
            # by absolute value of the lat coef
            col = sorted(
                enumerate(lat.column(my)),
                key=lambda mx_v: abs(mx_v[1]),
                reverse=True,
            )
            # keep only nonzero entries
            col = [mx_v for mx_v in col if mx_v[1]]
            columns.append(col)

        # cx defines which rows to flip
        #cnt_check = 0
        #n_cols = 0
        for cx in range(2**n):
            if not exclude_cx_cy:
                pairs = []
                for my in range(2**m):
                    # n_cols += 1
                    pos = neg = 0
                    for mx, value in columns[my]:
                        # cnt_check += 1
                        if parity_mx[mx & cx]:
                            value = -value
                        if not pos and value > 0:
                            pos = value
                        if not neg and value < 0:
                            neg = -value
                        if pos and neg:
                            break
                    pairs.append((pos, neg))

                for cy, val in enumerate(hadamard(pairs)):
                    if flip_input_output:
                        result[cy,cx] += val
                    else:
                        result[cx,cy] += val
            else:
                for cy in range(2**m):
                    val = 0
                    for my in range(2**m):
                        best = float("-inf")
                        for mx, value in columns[my]:
                            if parity_mx[mx & cx] ^ parity_my[my & cy]:
                                value = -value

                            #value <<= 1
                            if abs(value) + 1 <= best:
                                break

                            if flip_input_output:
                                if self.S[cy] == cx:
                                    value -= 1
                            else:
                                if self.S[cx] == cy:
                                    value -= 1
                            # if parity_mx[mx & cx] ^ parity_my[my ^ cy]:
                            #     value += 1
                            # else:
                            #     value -= 1

                            best = max(best, value)
                        val += best

                    if flip_input_output:
                        result[cy,cx] += val
                    else:
                        result[cx,cy] += val

        # print("cnt_check", cnt_check, "n_cols", n_cols, "avg", cnt_check / n_cols * 1.0)
        # ~3 for random 8-bit S-box
        return result


def hadamard(pairs, test=False):
    """For every u=0 to 2^n compute sum_x pairs[x][<u,x>].

    E.g. result[0] = sum(pair[0] for pair in pairs),
         result[1] chooses pair[0] for even indices and pair[1] for odd indices, etc.
    """
    n = int(len(pairs)).bit_length() - 1
    assert len(pairs) == 1 << n

    base, diff = [], []

    # convert into  base ± diff
    for a, b in pairs:
        # ensure integer average
        a <<= 1
        b <<= 1
        base.append((a+b)>>1) # average
        diff.append(a - base[-1]) # difference
        assert a == base[-1] + diff[-1]
        assert b == base[-1] - diff[-1]

    # Hadamard transform
    for i in range(n):
        step = 2**i
        for x in range(2**n):
            x2 = x ^ step
            if x < x2:
                a, b = diff[x], diff[x2]
                diff[x], diff[x2] = a + b, a - b

    s = sum(base)
    # remove scaling by 2
    ret =  [(s + d) >> 1 for d in diff]

    if test:
        test = []
        for cx in range(2**n):
            test.append(
                sum(
                    pair[Bin(cx, n).scalar_bin(x)]
                    for x, pair in enumerate(pairs)
                )
            )
        assert ret == test
    return ret


def test_hadamard():
    assert hadamard([(10, 20), (1, 2), (3, 4), (5, 6)], test=True) == [19, 21, 21, 21]
    assert hadamard([(10, 20), (1, 2), (31, 41), (52, 68)], test=True) == [94, 111, 120, 105]

    for n in range(8):
        print("hadamard n", n)
        for _ in range(10):
            lst = [
                (random.randrange(-1000, 1000), random.randrange(-1000, 1000))
                for _ in range(2**n)
            ]
            hadamard(lst, test=True)


def test_is_affine():
    # AES S-box
    S = (99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125, 250, 89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165, 229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117, 9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252, 177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127, 80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205, 12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144, 136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145, 149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186, 120, 37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14, 97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206, 85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22)
    alg = Algorithm(S, 8, 8)

    approx18 = (0, 1, 8, 29, 47, 53, 54, 57, 64, 74, 99, 102, 171, 179, 194, 211, 232, 239)
    assert alg.is_affine_approximation(approx18)
    assert alg.is_affine_approximation(approx18[:-1], check_maximal=False)
    assert not alg.is_affine_approximation(approx18[:-1], check_maximal=True)
    assert not alg.is_affine_approximation(approx18 + (2,))
    assert not alg.is_affine_approximation(approx18 + (2,), check_maximal=False)

    assert alg.expand_approximation(approx18) == approx18
    assert alg.expand_approximation(approx18[:17]) == approx18
    assert alg.expand_approximation(approx18[:14]) == approx18
    assert alg.expand_approximation(approx18[5:]) == approx18

    A = alg.approximation_from_inputs(approx18)

    # check it agrees with the approximation
    for x in approx18:
        assert A[x] == S[x]

    # check it's affine
    for _ in range(100):
        xs = random.sample(alg.domain, 3)
        xs.append(XOR(xs))
        assert XOR(A[x] for x in xs) == 0


def test_kzerosums():
    n_full_tests = 0
    n_itr = 0
    for n in range(3, 7):
        for m in range(1, 7):
            for itr in range(n_itr):
                print("itr", n, m, itr)
                if n == m and itr % 2:
                    S = list(range(2**n))
                    random.shuffle(S)
                else:
                    S = [random.randrange(2**m) for _ in range(2**n)]

                Alg = Algorithm(S, n, m)

                L4 = list(Alg.iter_kzerosums(4))
                for a, b, c, d in L4:
                    assert XOR([a,b,c,d]) == XOR([S[a],S[b],S[c],S[d]]) == 0

                L4_test = []
                for a, b, c, d in itertools.combinations(list(range(2**n)), 4):
                    if XOR([a,b,c,d]) == XOR([S[a],S[b],S[c],S[d]]) == 0:
                        L4_test.append((a,b,c,d))

                assert sorted(L4) == sorted(L4_test)

                for k in range(4, 11):
                    for _ in range(10):
                        if binomial(2**n, k) / QQ(2**(n+m)) <= 250:
                            continue
                        if 1 * 2**(k-1) > 2**n:
                            continue
                        print("test nmk", n, m, k)
                        xs = Alg.sample_kzerosum(k, affine_basis=True)
                        assert XOR(xs) == XOR(S[x] for x in xs) == 0
                        if k % 2 == 0:
                            assert len(AffineSpan(xs)) == 2**(k-2)
                        else:
                            assert len(AffineSpan(xs)) == 2**(k-1)
                        n_full_tests += 1
    print(n_full_tests)
    assert n_full_tests == 350 * n_itr


def test_lat_xors():
    n = 7
    S = list(range(2**n))
    random.shuffle(S)
    Alg = Algorithm(S, 7, 7)

    result1 = Alg._lat_xors_max_output()
    result2 = Alg.lat_xors(score_inputs=False)

    for score, cx, cy in result1:
        assert result2[cx,cy] == score + 2**(n-1)
