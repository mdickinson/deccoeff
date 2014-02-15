# Figure out best strategy for mixed-size Karatsuba.

# For an m-by-n multiplication:  choose 0 <= i <= m and 0 <= j <= n.
# Then Karatsuba allows us to split the multiplication into:
#    an i-by-j multiplication, and an
#    (m-i)-by-(n-j) multiplication, and an
#    max(i, m-i) by max(j, n-j) multiplication.

# ex: 2 by 2n -> 1 x n, 1 x n, 1 x n

# ex 2 by 20: compute 


# If i == 0 we get m-by-(n-j), plus m by max(j, n-j).  Better to
# do a simple split.  So include that possibility.

# Wlog assume 0 <= 2i <= m.

def minima(f, inputs):
    """Find minima of a function applied to various inputs.
    Return the inputs that minimize the function, and
    the minimum value."""

    min = None
    best = []
    for x in inputs:
        v = f(x)
        if min is None or v < min:
            min = v
            best = [x]
        elif v == min:
            best.append(x)
    return min, best

def Kmemo(x, _memo = {}):
    if x not in _memo:
        _memo[x] = K(x)
    return _memo[x]


def K0(m, n):
    return K(m, n)[0]

def K(m, n, _memo = {}):
    """Given m and n, return the minimum number of multiplications
    needed, and the splits that achieve that minimum."""
    assert(m >= 1 and n >= 1)

    if (m, n) in _memo:
        return _memo[m, n]

    if m == n == 1:
        _memo[m, n] = 1, []
    else:
        def f(i, j):
            assert not ((i == 0 or i == m) and (j == 0 or j == n))
            if i == 0 or i == m:
                return K0(m, j) + K0(m, n-j)
            elif j == 0 or j == n:
                return K0(i, n) + K0(m-i, n)
            else:
                return (K0(i, j) + K0(m-i, n-j) +
                        K0(max(i, j, m-i), max(i, j, n-j)))

        splits = [(i, j)
                  for i in range(m+1)
                  for j in range(n+1)
                  if not ((i == 0 or i == m) and (j == 0 or j == n))]

        _memo[m, n] = minima(lambda x: f(x[0], x[1]), splits)
    return _memo[m, n]

# Here's one possible optimal(?) strategy

def next_higher_power_of_two(m):
    i = 1
    while m > i:
        i *= 2
    return i

def opt(m, n):
    """Return number of multiplications required."""
    assert(m >= 1 and n >= 1)
    if n < m:
        return opt(n, m)

    if m == n == 1:
        # base case: single multiplication
        return 1

    k = next_higher_power_of_two(m)
    if n > k:
        # straight split, m * n = m * k + m * (n-k)
        # ex: 9 * 23 -> 9 * 16 + 9 * 7
        return opt(m, k) + opt(m, n-k)
    else:
        # skew Karatsuba:
        #  split 1st arg as [0:k/2], [k/2:m]
        #  split 2nd arg as [0:n-k/2], [n-k/2:n]
        return opt(k//2, n-k//2) + \
            opt(m-k//2, k//2) + \
            opt(k//2, k//2)

for m in range(1, 100):
    for n in range(1, 100):
        print m, n
        K(m, n)
