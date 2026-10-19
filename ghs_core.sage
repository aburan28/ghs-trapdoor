# Core field utilities for the GHS trapdoor demo (n=4, l=5: K=F_{2^20}, k=F_{2^5})
l = 5; n = 4; N = l*n
F2 = GF(2)
K.<a> = GF(2^N, 'a')
k = K.subfield(l, 'bb'); b = k.gen()
g = K.multiplicative_generator()
theta = g^((2^N-1)//5)          # primitive 5th root of unity, Tr_{F16/F2}(theta)=1
th = [theta^(2^i) for i in range(4)]  # normal basis over k, sum = 1

def f2vec(x):
    v = list(x.polynomial()); return vector(F2, v + [0]*(N-len(v)))
def f2_from_vec(v):
    return K({i: v[i] for i in range(N) if v[i]})

# F2-basis {b^s * th[j]}
_base = [b^s * th[j] for j in range(4) for s in range(l)]
_B = Matrix(F2, [f2vec(e) for e in _base]).transpose()
assert _B.rank() == N
_Binv = _B.inverse()

def normal_coords(x):
    """K-element -> [y_0..y_3] with x = sum y_j*th[j], y_j in k"""
    c = _Binv * f2vec(x)
    return [sum(c[l*j+s]*b^s for s in range(l)) for j in range(4)]

def from_normal(ys):
    return sum(ys[j]*th[j] for j in range(4))

def sigma(x):
    return x^(2^l)

def sqrtK(x):
    return x^(2^(N-1))
