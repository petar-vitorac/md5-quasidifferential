import pyboolector
import sys
import math
import numpy as np

##################
# INITIALIZATION #
##################

# Increase recursion limit to be able to compute REF of big matrix
sys.setrecursionlimit(3000)

BITS_PER_WORD = 32

# Only use the first n steps (entire md5 = 64 steps)
nb_steps = 64

# Initialize Boolector
btor = pyboolector.Boolector()
btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)
btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)

# Initialize bit vector variables for masks on registers (and branches)
uA = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uA%d" % i) for i in range(1, nb_steps+2)]
uB = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uB%d" % i) for i in range(1, nb_steps+2)]
uC = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uC%d" % i) for i in range(1, nb_steps+2)]
uD = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uD%d" % i) for i in range(1, nb_steps+2)]
uB1 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uB1_%d" % i) for i in range(1, nb_steps+1)]
uB2 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uB2_%d" % i) for i in range(1, nb_steps+1)]
uC1 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uC1_%d" % i) for i in range(1, nb_steps+1)]
uD1 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uD1_%d" % i) for i in range(1, nb_steps+1)]

# The masks on the state after the last step are zero
btor.Assert(uA[nb_steps+1] == 0)
btor.Assert(uB[nb_steps+1] == 0)
btor.Assert(uC[nb_steps+1] == 0)
btor.Assert(uD[nb_steps+1] == 0)

# Wang's characteristic
char = [[]]*8 + [list(range(6,23)), [6,23,31], list(range(0,12))+list(range(23,32)), [0,23]+list(range(15,21)), [0,1,6,7,8,31], [12,13,31], [30,31], [7,8,31]+list(range(13,20)), [24,25,31], [31], [3,15,31], [29,31]] + [[31]]*2 + [[17,31]] + [[31]]*3 + [[]]*12 + [[31]]*27 + [[25,31], [25,26,31], [25,31]]

# Initialize sign expression
sign = btor.Const(0, 1)

# The differences on the registers
aBs = [btor.Const(sum(map(lambda x: 2**x, el)), BITS_PER_WORD) for el in char]
#aBs = aBs[:len(uBs)]
# aX(i) is the difference on register X at the i-th step (numbering from 1)
aB = lambda i : aBs[i+2]
aC = lambda i : aBs[i+1]
aD = lambda i : aBs[i]
aA = lambda i : aBs[i-1]

# The constatnts K_i and masks on them
K = [None] + [btor.Const(int(2**32*abs(math.sin(i+1))), BITS_PER_WORD) for i in range(nb_steps)]
uK = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uK_%d" % i) for i in range(1, nb_steps+1)]

# The number of bit shifts after the third addition
s = [None] + [btor.Const(i, BITS_PER_WORD) for i in [7, 12, 17, 22]*4 + [5, 9, 14, 20]*4 + [4, 11, 16, 23]*4 + [6, 10, 15, 21]*4]

# The differences and masks on the message
aM = [None] + [btor.Const(i, BITS_PER_WORD) for i in [0]*4 + [2**31] + [0]*6 + [2**15] + [0]*2 + [2**31] + [0]]
uM = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "uM_%d" % i) for i in range(1, nb_steps+1)]

# UNCOMMENT THE CHARACTERISTIC YOU WISH TO USE #
# Intermediary differences for collision from Wang's paper
b1 = [0, 0, 0, 0, 0, 526336, 8371200, 167838756, 2164327744, 2224029761, 8396865, 257, 2148401280, 2148524032, 2147745792, 2181038080, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 0, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 0, 2147483648, 2147483648, 2147483648]
b2 = [0, 0, 0, 0, 0, 526336, 115712, 1979909220, 2264990976, 202375171, 402661376, 8618752, 131263, 4096, 1077673984, 101179776, 16777216, 0, 32776, 536870912, 0, 0, 917504, 2147483648, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 0, 0]
b3 = [0, 0, 0, 0, 2147483648, 7878656, 117760, 167976420, 2198277888, 202375169, 135258112, 8782080, 131151, 28672, 3221487616, 503341184, 117440512, 0, 8, 536870912, 0, 0, 393216, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32768, 0, 0]
b4 = [0, 0, 0, 0, 2147483648, 526336, 48128, 1979780068, 2164718848, 68157443, 134225920, 33423616, 131139, 4096, 1074528256, 234905728, 16777216, 0, 120, 3758096384, 0, 0, 131072, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32768, 0, 0]

# Intermediary differences for collision from computed collision 1
"""
b1 = [0, 0, 0, 0, 0, 526336, 8371200, 167838756, 2164327744, 2224029761, 8396865, 257, 2148401280, 2148524032, 2147745792, 2181038080, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 0, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 2147483648, 0, 2147483648, 2147483648, 2147483648]
b2 = [0, 0, 0, 0, 0, 526336, 52224, 437193764, 2264990976, 202375175, 134225920, 25395968, 131139, 4096, 3225157632, 33564544, 16777216, 0, 98312, 536870912, 0, 0, 131072, 2147483648, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2147483648, 0, 0, 0]
b3 = [0, 0, 0, 0, 2147483648, 7878656, 248832, 167976420, 2198277888, 68157441, 134471680, 8785664, 131265, 126976, 1074528256, 100688000, 251658240, 0, 8, 1610612736, 0, 0, 393216, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 98304, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 491520, 0, 0]
b4 = [0, 0, 0, 0, 2147483648, 526336, 48128, 1979780068, 2164391168, 68157443, 134225920, 27132160, 131137, 4096, 1075576832, 33562752, 50331648, 0, 24, 536870912, 0, 0, 131072, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32768, 0, 0]
"""

# masks and differences at the output of the F function
v1 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "v1_%d" % i) for i in range(1, nb_steps+1)]
b1 = [None] + [btor.Const(i, BITS_PER_WORD) for i in b1]

# masks and differences at the output of the first modular addition (the one with F)
v2 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "v2_%d" % i) for i in range(1, nb_steps+1)]
b2 = [None] + [btor.Const(i, BITS_PER_WORD) for i in b2]

# masks and differences at the output of the second modular addition (the one with the message)
v3 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "v3_%d" % i) for i in range(1, nb_steps+1)]
b3 = [None] + [btor.Const(i, BITS_PER_WORD) for i in b3]

# masks and differences at the output of the third modular addition (the one with the K_i constants)
v4 = [None] + [btor.Var(btor.BitVecSort(BITS_PER_WORD), "v4_%d" % i) for i in range(1, nb_steps+1)]
b4 = [None] + [btor.Const(i, BITS_PER_WORD) for i in b4]

##########################
# BEGIN HELPER FUNCTIONS #
##########################

# M pseudoinverse, as used for computing the QDTM for modular addition
def M_pseudoinverse(t):
    t = t ^ (t << 1)
    return (t >> 1)

# M transpose, as used for computing the QDTM for modular addition
def M_transpose(t):
    t = (t >> 1)
    i = 1
    while i < BITS_PER_WORD:
        t = t ^ (t >> i)
        i *= 2
    return t

# Hamming weight of a 32-bit vector
def wt(x):
    x -= (x >> 1) & 0x55555555
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0F0F0F0F
    x += (x >> 8)
    x += (x >> 16)
    return x & 0x0000003F

# Parity of the dot product of vectors
def dot(a, b):
    return wt(a&b)[0]

# Assert what is necessary for going through the f function, update the sign, and return the weight
def f_asserts(btor, a, b, c, d, u, v, w, t):
    global sign
    sign ^= dot((a&b) ^ (~a&~b&u), v) ^ dot((a&v) ^ (~a&b&u), c^d)

    btor.Assert(t == w)
    btor.Assert((a & u) ^ (b & v) == t & (c ^ d))
    btor.Assert((c^d) & ~(a|b) == 0)
    btor.Assert((u|v) & ~(t|a|b) == 0)

    return wt(t | a | b)


# Assert what is necessary for going through a modular addition, update the sign, and return the weight (partly from https://github.com/TimBeyne/quasidifferential-trails)
def add_asserts(btor, a, b, c, u, v, w):
    global sign
    nb_bits = BITS_PER_WORD

    a_ = b ^ c
    b_ = a ^ c
    c_ = M_pseudoinverse(a ^ b ^ c)

    u_ = u ^ w
    v_ = v ^ w
    w_ = M_transpose(u ^ v ^ w)

    n = nb_bits - 1

    sign ^= dot((~a_ & u_) ^ (c_ & v_), (~b_ & v_) ^ (c_ & u_)) ^ dot(u_ & v_, c_ ^ (a_ & b_ & ~c_))

    btor.Assert((u_ | v_) & ~(a_ | b_ | w_) == 0)
    btor.Assert((a_ & u_) ^ (b_ & v_) == (c_ & w_))
    btor.Assert(((a_[n] == 0) & (b_[n] == 0)) | (a_[n] & u_[n] == u_[n] ^ v_[n]))

    weight = wt(w_ & ~a_ & ~b_)
    extra = btor.Cond(
        (a_[n] | b_[n]) & ((u_[n] ^ v_[n]) == (a_[n] & u_[n])) & (u_[n] != v_[n]),
        btor.Const(1, nb_bits),
        btor.Const(0, nb_bits)
    )
    return weight - extra



# For finding the naively computed average weight of the characteristic
def assert_masks_0():
    for x in uA[1:]:
        btor.Assert(x == 0)
    for x in uB[1:]:
        btor.Assert(x == 0)
    for x in uB1[1:]:
        btor.Assert(x == 0)
    for x in uB2[1:]:
        btor.Assert(x == 0)
    for x in uC[1:]:
        btor.Assert(x == 0)
    for x in uC1[1:]:
        btor.Assert(x == 0)
    for x in uD[1:]:
        btor.Assert(x == 0)
    for x in uD1[1:]:
        btor.Assert(x == 0)
    for uKi in uK[1:]:
        btor.Assert(uKi == 0)
    for uMi in uM[1:]:
        btor.Assert(uMi == 0)
    for v1i in v1[1:]:
        btor.Assert(v1i == 0)
    for v2i in v2[1:]:
        btor.Assert(v2i == 0)
    for v3i in v3[1:]:
        btor.Assert(v3i == 0)
    for v4i in v4[1:]:
        btor.Assert(v4i == 0)
# assert_masks_0() # Uncomment this to only find the quasidifferential trail with null mask

# Convert masks on messages in individual steps to the global message mask
def get_uM_final(u):
    final = list()
    for i in range(min(16, nb_steps)):
        final.append(u[i])
    for i in range(16, min(32, nb_steps)):
        final[(5*i+1)%16] ^= u[i]
    for i in range(32, min(48, nb_steps)):
        final[(3*i+5)%16] ^= u[i]
    for i in range(48, min(64, nb_steps)):
        final[(7*i)%16] ^= u[i]
    return final

# Helper function for binary manipulation
def bin_idx(n, idx):
    return int(bin(n)[-idx-1])

# Helper function for binary manipulation
def bin_size(n):
    return len(bin(n))-2

# Helper function that returns the non null bits in a mask
def non_null_bits(u):
    non_null = list()
    for ui in u:
        non_null.append(list())
        for i in range(bin_size(ui)):
            if bin_idx(ui, i):
                non_null[-1].append(i)
    return non_null

# Assert that the future solution is not equal to one of the previously computed solutions
def distinct_assumptions(btor, masks, solutions):
    result = btor.Const(1)
    for solution in solutions:
        here = btor.Const(0)
        for i in range(len(masks)):
            for j in range(len(masks[i])):
                here |= (masks[i][j] != btor.Const(solution[i][j], BITS_PER_WORD))
        result &= here
    btor.Assume(result)

# Row echelon form of a binary matrix
def ref(A):
    r, c = A.shape
    if r == 0 or c == 0:
        return A

    for i in range(len(A)):
        if A[i, 0] != 0:
            break
    else:
        B = ref(A[:,1:])
        return np.hstack([A[:,:1], B])
    
    if i > 0:
        ith_row = A[i].copy()
        A[i] = A[0]
        A[0] = ith_row

    for i in range(1, r):
        if A[i, 0]:
            A[i] = np.logical_xor(A[i], A[0])

    B  = ref(A[1:, 1:])

    return np.vstack([A[:1], np.hstack([A[1:, :1], B])])

# Assert that the future solution is not linear combination of previously computed solutions
def not_linear_comb_assumptions(btor, this, solutions):
    if not len(solutions):
        return
    solutions = [list(map(int, list(''.join(sum(solution, []))))) for solution in solutions]
    A = np.array(solutions).transpose()
    B = np.hstack([A, np.eye(A.shape[0])])
    CT = ref(B)
    C = CT[:,:A.shape[1]]
    T = CT[:,A.shape[1]:]
    this = sum(this, [])

    found_one = btor.Const(0, 1)
    for i in range(len(solutions), B.shape[0]):
        expr = btor.Const(0, 1)
        # take the row from T
        Ti = T[i,:]
        
        for j, mask in enumerate(this):
            Ti_local = Ti[j*BITS_PER_WORD:(j+1)*BITS_PER_WORD]
            Ti_local = Ti_local.dot(1 << np.arange(Ti_local.size)[::-1])
            expr ^= dot(mask, btor.Const(int(Ti_local), BITS_PER_WORD))

        found_one |= expr
    
    btor.Assume(found_one == btor.Const(1, 1))

##################################
# BEGIN MAIN PART OF THE PROGRAM #
##################################

# Going through each step, adding all the assertions and calculating the weight of the trail
weight = btor.Const(0, BITS_PER_WORD)
for i in range(1, nb_steps+1):
    # Asserts for the branching
    btor.Assert(uB[i] == uB1[i] ^ uB2[i] ^ uC[i+1])
    btor.Assert(uC[i] == uC1[i] ^ uD[i+1])
    btor.Assert(uD[i] == uD1[i] ^ uA[i+1])

    # Asserts for the F functions
    if i <= 16:
        weight += f_asserts(btor,
            aB(i), (aC(i) ^ aD(i)), aD(i), b1[i],
            uB1[i], uC1[i], (uC1[i] ^ uD1[i]), v1[i])
        g = i
    elif i <= 32:
        weight += f_asserts(btor,
            (aB(i) ^ aC(i)), aD(i), aC(i), b1[i],
            uB1[i], uD1[i], (uB1[i] ^ uC1[i]), v1[i])
        g = (5*(i-1)+1)%16 + 1
    elif i <= 48:
        btor.Assert(b1[i] == aB(i) ^ aC(i) ^ aD(i))
        btor.Assert(v1[i] == uB1[i] == uC1[i] == uD1[i])
        g = (3*(i-1)+5)%16 + 1
    else:
        weight += f_asserts(btor,
            aB(i), aD(i), (aC(i) ^ aD(i)), b1[i],
            uB1[i], (uC1[i] ^ uD1[i]), uC1[i], v1[i])
        sign ^= wt(v1[i])[0]
        g = (7*(i-1))%16 + 1

    # Asserts for the modular addition
    weight += add_asserts(btor, b1[i], aA(i), b2[i], v1[i], uA[i], v2[i])
    weight += add_asserts(btor, b2[i], aM[g], b3[i], v2[i], uM[i], v3[i])
    weight += add_asserts(btor, b3[i], btor.Const(0, BITS_PER_WORD), b4[i], v3[i], uK[i], v4[i])
    weight += add_asserts(btor, aB(i), btor.Rol(b4[i], s[i]), aB(i+1), uB2[i], btor.Rol(v4[i], s[i]), uB[i+1])


# Calculate "true sign" (i.e. remove dependence on the mask on constants K. In other words, dot(uM, message) = (-1)^true_sign will be the condition)
true_sign = sign
for i in range(1, nb_steps+1):
    true_sign ^= dot(uK[i], K[i])

# Calculate the global mask on the message (combining the per-step masks)
uM_final = get_uM_final(uM[1:])

# Parameters for trail searching
min_target_weight = 0 # start by searching trails with correlation 2^(-min_target_weight)
max_target_weight = 222 # stop searching trails when reached correlation 2^(-max_target_weight)
max_solutions_per_weight = 50 # don't show more than max_solutions_per_weight trails per weight

solutions = []
sign_solutions = []

# Start at min_target_weight, and progressively increment weight to find most important trails first
for target_weight in range(min_target_weight, max_target_weight+1):
    print("trying weight " + str(target_weight))
    current_solutions = 0
    while current_solutions < max_solutions_per_weight:
        btor.Assume(weight == target_weight)
        distinct = [uM_final]
        # distinct_assumptions(btor, distinct, solutions)
        not_linear_comb_assumptions(btor, distinct, solutions)
        to_print = [uM[1:], uK[1:], uM_final]
        r = btor.Sat()
        if r == btor.SAT:
            sol = [[xi.assignment for xi in x] for x in distinct]
            solutions.append(sol)
            sign_solutions.append(true_sign.assignment)
            current_solutions += 1
            print("Mask on the message (. represents 0 mask):")
            for xi in uM_final:
                if int(xi.assignment, base = 2) == 0:
                    print('.', end='')
                else:
                    print(hex(int(xi.assignment, base = 2)), end=' ')
            print()
            print(f"\"True\" sign: {true_sign.assignment}")
            print("Non null bits on message mask:")
            print(non_null_bits([int(uMi.assignment, base=2) for uMi in uM_final]))
            print()
            print(f"Found {current_solutions} solutions so far")
        else:
            print(f"{current_solutions} solutions found")
            break

    if current_solutions == max_solutions_per_weight:
        print(f"Stopped because {current_solutions} found for this weight")

