from sage.crypto.sbox import SBox
import functools
import itertools

# S-box of the example cipher
S = SBox([7, 2, 4, 5, 1, 6, 3, 0])

# Bit permutations used in the example cipher
def permute(x):
    y = vector(GF(2), 9)
    for i in range(9):
        s = i // 3 # S-box index
        j = i % 3  # Bit

        y[8 - (3 * j + ((s + 1) % 3))] = x[8 - i]
    return y

# Helper function for indexing in the quasidifferential transition matrix
def to_index(mask):
    return ZZ(list(mask)[::-1], base=2)

# Helper function for calculating the quasidifferential transition matrix of the S-Box (https://github.com/TimBeyne/quasidifferential-trails)
def interleave_bits(x, y, n):
    z = 0
    for i in range(n):
        z |= (x & (1 << i)) << i | (y & (1 << i)) << (i + 1);
    return z

# Helper function for calculating the quasidifferential transition matrix of the S-Box (https://github.com/TimBeyne/quasidifferential-trails)
def to_quasidifferential_basis(x):
    if len(x) == 1:
        return x

    assert len(x) % 4 == 0

    l = int(len(x) / 4) 

    x_00 = to_quasidifferential_basis(x[   :  l])
    x_01 = to_quasidifferential_basis(x[  l:2*l])
    x_10 = to_quasidifferential_basis(x[2*l:3*l])
    x_11 = to_quasidifferential_basis(x[3*l:   ])

    return vector(x.base_ring(),
        (x_00 + x_11).list() +\
        (x_01 + x_10).list() +\
        (x_00 - x_11).list() +\
        (x_01 - x_10).list())

# Helper function for calculating the quasidifferential transition matrix of the S-Box (https://github.com/TimBeyne/quasidifferential-trails)
def interleaved_transition_matrix(F, n, m):
    T = matrix(QQ, 2 ** (2*m), 2 ** (2*n))
    for x in range(2 ** n):
        for y in range(2 ** n):
            i = interleave_bits(x, y, n)
            j = interleave_bits(F(x), F(y), m)
            T[j, i] = 1
    return T

# Calculate the quasidifferential transition matrix of the S-Box (step 1) (https://github.com/TimBeyne/quasidifferential-trails)
def quasidifferential_transition_matrix(F, n, m):
    D = interleaved_transition_matrix(F, n, m)

    for i in range(2 ** (2*n)):
        D.set_column(i, to_quasidifferential_basis(D.column(i)))
    for i in range(2 ** (2*m)):
        D.set_row(i, to_quasidifferential_basis(D.row(i)))

    return D / 2**n

# Calculate the quasidifferential transition matrix of the S-Box (step 2) (https://github.com/TimBeyne/quasidifferential-trails)
def deinterleave_quasidifferential_transition_matrix(D, n, m, primary = 'diff'):
    R = matrix(QQ, 2 ** (2*m), 2 ** (2*n))
    for u, v in itertools.product(range(2 ** n), range(2 ** m)):
        for a, b in itertools.product(range(2 ** n), range(2 ** m)):
            if primary == 'mask':
                R[2 ** m * v + b, 2 ** n * u + a] = D[
                    interleave_bits(b, v, m), interleave_bits(a, u, n)
                ] 
            elif primary == 'diff':
                R[2 ** m * b + v, 2 ** n * a + u] = D[
                    interleave_bits(b, v, m), interleave_bits(a, u, n)
                ] 
    return R

# Call the functions and store the QDTM in DF
DF = quasidifferential_transition_matrix(S, 3, 3)
DF = deinterleave_quasidifferential_transition_matrix(DF, 3, 3)


"""
"Bruteforce" way to calculate QDT matrix:
For each input/output mask, iterate over all right pairs and calculate how many satisfy the relation.

part_DFs = [[None for x in range(2**3)] for y in range(2**3)] 
def part_DF(a, b):
    if part_DFs[to_index(a)][to_index(b)]:
        return part_DFs[to_index(a)][to_index(b)]
    
    xs = [x for x in VectorSpace(GF(2), 3) if S(x+a) == S(x) + b]
    result = matrix(QQ, 8, 8, 0)
    if len(xs) == 0:
        return result

    for v in VectorSpace(GF(2), 3):
        for u in VectorSpace(GF(2), 3):
            count = len([x for x in xs if u.dot_product(x) == v.dot_product(S(x))])
            result[to_index(v), to_index(u)] = len(xs)/8 * (2*count/len(xs)-1)

    part_DFs[to_index(a)][to_index(b)] = result
    return result
"""

# Return a function that only queries the QDTM for a certain input and output difference a and b.
def part_DF(a, b):
    def f(u, v):
        return DF[to_index(a)*2**3+to_index(u), to_index(b)*2**3+to_index(v)]
    
    return f
    
# Characteristic
diff_trails = [
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 0, 1]), vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])],
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 1, 0, 1]), vector(GF(2), [0, 1, 0, 0, 0, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])],
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 0, 1, 1]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])],
    [vector(GF(2), [0, 0, 1, 0, 0, 0, 0, 0, 0]), vector(GF(2), [0, 0, 0, 0, 0, 0, 1, 1, 1]), vector(GF(2), [0, 1, 0, 0, 1, 0, 0, 1, 0]), vector(GF(2), [0, 0, 0, 0, 1, 0, 0, 0, 0])]
]

#roundkeys  = [vector(GF(2), 9) for _ in range(4)]
roundkeys  = [vector(GF(2), [1, 0, 1, 1, 0, 0, 0, 1, 1]) for _ in range(3)] + [vector(GF(2), 9)]

initial_mask = vector(GF(2), 9)
initial_mask.set_immutable()

curr_trail = list()

# Compute quasidifferential trails and their probability with "path-searching" algorithm
@functools.cache
def trail(d, pos, mask):
    global curr_trail
    differences = diff_trails[d]
    if len(curr_trail) == len(differences) - 1:
        return 1
    
    prev_mask = permute(mask)
    prev_diff = permute(differences[pos-1])
        
    masks = VectorSpace(GF(2), 9)
    if pos == len(differences)-1:
        masks = [initial_mask]
    
    result = 0
    
    S1_DFi = part_DF(differences[pos][0:3], prev_diff[0:3])
    S2_DFi = part_DF(differences[pos][3:6], prev_diff[3:6])
    S3_DFi = part_DF(differences[pos][6:9], prev_diff[6:9])
        
    for m in masks:
        DFi_here = S1_DFi(m[0:3], prev_mask[0:3]) * S2_DFi(m[3:6], prev_mask[3:6]) * S3_DFi(m[6:9], prev_mask[6:9])
        if DFi_here == 0:
            continue
        DFi_here *= (-1)**prev_mask.dot_product(roundkeys[pos-1])
        m.set_immutable()
        curr_trail += [m]
        if len(curr_trail) == len(differences) - 1:
            print(curr_trail)
        result += DFi_here * trail(d, pos+1, m)

        curr_trail.pop()
        
    return result

# Output should match that of cipher.sage
results = list()
for i in range(len(diff_trails)):
    results.append(trail(i, 1, initial_mask))
    print(results[-1])
print("Sum:", sum(results))

