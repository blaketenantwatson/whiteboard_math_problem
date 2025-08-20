from math import isqrt
from typing import Optional, Tuple, List, Dict
from tqdm import tqdm

def is_square(n: int) -> bool:
    if n < 0:
        return False
    s = isqrt(n)
    return s * s == n

def legendre_symbol(a: int, p: int) -> int:
    """
    Legendre symbol (a/p) mod odd prime p:
      returns 1 if a is a nonzero quadratic residue mod p
              0 if a ≡ 0 (mod p)
             -1 if a is a non-residue
    Implemented via Euler's criterion.
    """
    a %= p
    if a == 0:
        return 0
    # pow returns in [0, p-1]; map p-1 to -1 for convenience
    r = pow(a, (p - 1) // 2, p)
    return 1 if r == 1 else (-1 if r == p - 1 else 0)

def miller_rabin(n: int, bases: Optional[List[int]] = None) -> bool:
    """Deterministic MR for 64-bit, probabilistic otherwise (kept for single-k use)."""
    if n < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
    for p in small_primes:
        if n % p == 0:
            return n == p
    # write n-1 = d*2^s
    d = n - 1
    s = 0
    while (d & 1) == 0:
        d >>= 1
        s += 1
    if bases is None:
        if n < 2_152_302_898_747:
            bases = [2, 3, 5, 7, 11]
        elif n < 3_474_749_660_383:
            bases = [2, 3, 5, 7, 11, 13]
        elif n < 341_550_071_728_321:
            bases = [2, 3, 5, 7, 11, 13, 17]
        else:
            bases = [2, 3, 5, 7, 11, 13, 17, 19, 23]

    def check(a: int) -> bool:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            return True
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                return True
        return False

    for a in bases:
        if a % n and not check(a):
            return False
    return True

def representable_by_form(k: int, N: int, want_solution: bool = True) -> Optional[Tuple[int, int]]:
    """
    Check/solve x^2 + x*y + 2k*y^2 = N over integers.
    Returns (x,y) or None. O(sqrt(N/p)) loop where p = 8k-1.
    """
    if N < 0:
        return None
    p = 8 * k - 1
    if p <= 0:
        return None

    # y=0 -> x^2=N (fast path)
    if is_square(N):
        return (isqrt(N), 0) if want_solution else (0, 0)

    # Necessary condition mod p: 4N must be a square => (N/p)=1
    # Quick exit to avoid a pointless y-scan:
    if legendre_symbol(N, p) != 1:
        return None

    # Discriminant D = -p*y^2 + 4N must be a perfect square; bound: D>=0 => |y|<=sqrt(4N/p)
    fourN = N << 2  # 4N
    Ymax_denom = p
    if Ymax_denom <= 0:
        return None
    Ymax = isqrt(fourN // Ymax_denom)

    # full symmetric scan (safe & simple)
    for y in range(-Ymax, Ymax + 1):
        yy = y * y
        D = fourN - p * yy
        if D < 0:
            continue
        s = isqrt(D)
        if s * s != D:
            continue
        # parity: (-y + s) must be even for integer x
        if ((-y + s) & 1) != 0:
            continue
        x1 = (-y + s) >> 1
        if x1 * x1 + x1 * y + 2 * k * yy == N:
            return (x1, y) if want_solution else (0, 0)
        x2 = (-y - s) >> 1
        if x2 * x2 + x2 * y + 2 * k * yy == N:
            return (x2, y) if want_solution else (0, 0)
    return None

def find_m_for_k(k: int, m_limit: int = 10_000, require_p_prime: bool = True):
    """
    Find m s.t. form represents m^3 but NOT m.  
    Uses Legendre prefilter to prune ~1/2 of m upfront.
    """
    p = 8 * k - 1
    if require_p_prime and not miller_rabin(p):
        return None

    for m in tqdm(range(2, m_limit + 1)):
        # skip trivial y=0 representables
        if is_square(m):
            continue

        # Necessary condition for BOTH m and m^3:
        # (m/p) must be 1, otherwise neither can be represented.
        if legendre_symbol(m, p) != 1:
            continue

        # Now m might be representable; we need it to be NOT.
        if representable_by_form(k, m, want_solution=True) is not None:
            continue

        # Check m^3
        N3 = m * m * m
        if representable_by_form(k, N3, want_solution=True) is not None:
            return (m, None, None)
    return None

def sieve_primes_upto(limit: int) -> List[bool]:
    """
    Exact sieve of Eratosthenes: returns is_prime array of length limit+1.
    """
    if limit < 2:
        return [False] * (limit + 1)
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    # even multiples
    for x in range(4, limit + 1, 2):
        is_prime[x] = False
    r = isqrt(limit)
    for i in range(3, r + 1, 2):
        if is_prime[i]:
            step = i << 1
            start = i * i
            is_prime[start:limit + 1:step] = [False] * (((limit - start) // step) + 1)
    return is_prime

def generate_m_map(N: int, m_limit: int = 200, require_p_prime: bool = True) -> Dict[int, int]:
    """
    Return a dict mapping k -> m for 1<=k<=N. Uses sieve to
    select only k with p=8k-1 prime (exact, fast for large N).
    """
    out: Dict[int, int] = {}
    if require_p_prime:
        max_p = 8 * N - 1
        is_prime = sieve_primes_upto(max_p)  # exact primality for all p in one pass
        ks = [k for k in range(1, N + 1) if is_prime[8 * k - 1]]
    else:
        ks = list(range(1, N + 1))

    for k in ks:
        res = find_m_for_k(k, m_limit=m_limit, require_p_prime=False)
        if res is not None:
            out[k] = res[0]
    return out

def sanity_check_map(m_map: Dict[int, int]) -> None:
    """
    Wasn't trusting the outputs, so this is quick sanity check on the outputs.
    Only works the map, you gotta check it manually for single values, or just wrap it in a map.
    Only prints when a sanity check fails for some (k,m).
    """
    for k, m in sorted(m_map.items()):
        has_m = representable_by_form(k, m, want_solution=True)
        has_m3 = representable_by_form(k, m**3, want_solution=True)
        if has_m is not None or has_m3 is None:
            print(f"[FAIL] k={k}, m={m}: has_m={'YES' if has_m else 'NO'}, has_m3={'YES' if has_m3 else 'NO'}")
    print("Finished sanity checks.")

if __name__ == "__main__":
    # Example map generation and sanity check
    m_map = generate_m_map(10000, m_limit=5000)
    for k, m in m_map.items():
        print(f"k = {k}, m = {m}")

    # Sanity check only prints problems
    sanity_check_map(m_map)

    # Example usage of single m value calculation
    k = 16
    result = find_m_for_k(k, m_limit=1000)
    print(result[0])
