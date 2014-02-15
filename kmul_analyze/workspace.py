# Analyze workspace required for Karatsuba algorithm

# Here's the skeleton of the algorithm, that does everything except
# the actual arithmetic.

def kmul(m, n, k):
    """Karatsuba multiplication, m-by-n.  Returns number
    of multiplications, and workspace required."""
    assert k < m <= 2*k and k < n <= 2*k and is_power_of_two(k)
    m1, w1 = dmul(k, k)   # workspace also needed for size-2k result
    m2, w2 = dmul(m-k, k) # result goes straight in..
    m3, w3 = dmul(n-k, k) # ditto
    return sum(m1, m2, m3), 2*k + max(w1, w2, w3)

def dmul(m, n):
    """Dispatch multiplication."""
    assert 1 <= m <= n
    if m == 1:
        # basecase multiplication; n multiplications, no workspace required
        return n, 0

    # twok = next power of 2 after m (possibly equal to m)
    twok = 1
    while twok < m:
        twok *= 2
    assert twok >= 2
    k = twok/2

    if n <= twok:
        # balanced: do single Karatsuba
        return kmul(m, n, k)

    # otherwise, split off an m * twok multiplication
    m1, w1 = kmul(m, twok, k)
    n -= twok

    # w = workspace required thus far
    w = w1

    while (n > twok):
        # need workspace to do this multiplication, *and* store the result
        m2, w2 = kmul(m, twok, k)
        w = max(w, w2 + m + twok)
        n -= twok

    # left with m-by-n multiplication; again need workspace + result
    m3, w3 = dmul(m, n)
    w = max(w, w3 + m + n)

    return sum(m1, m2, m3), w
