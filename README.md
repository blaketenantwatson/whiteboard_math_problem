# whiteboard_problem

# Quadratic Form Solver and Optimizer

Take everything in this with a grain of salt, this is just a first attempt at doing this. I'm no math major, I just am messing with this randomly. Most of these optimizations and math concepts are ripped from some combination of math overflow or LLMs.

I'm using `tqdm` to monitor progress, but if you don't feel like `pip` installing it, just remove the wrapper in `find_m_for_k` and the import line.

This project finds integers `m` such that for a given `k` (with `p = 8k - 1` prime):

- `x^2 + x*y + 2k*y^2 = m^3` **has** integer solutions, and  
- `x^2 + x*y + 2k*y^2 = m` **does not**.

---

## General Approach

1. **Mathematical Basis**  
   The equation `x^2 + x*y + 2k*y^2` is the *norm form* of the quadratic field  
   `Q(√(-p))`. Whether `N` is representable depends on whether it is a norm in this field.

2. **Algorithm Outline**  
   - For each candidate `m`, check:
     - Is `m` representable? (skip if yes)  
     - Is `m^3` representable? (accept if yes)  
   - Return the first such `m` (or `None` if none found within limits).

---

## Speed Optimizations

1. **Legendre Symbol Prefilter**  
   - If `(N/p) ≠ 1`, then `N` cannot be represented.  
   - This eliminates ~50% of candidates early.

2. **Residue Class Sieving for `y`**  
   - Instead of scanning every `y` up to `sqrt(4N/p)`, only scan congruence classes `y mod M0` that could possibly make the discriminant a square.  
   - Using moduli like `[8, 3, 5, 7, 11, 13]` prunes ~95% of impossible `y`.

3. **Prime Sieving for Batch Mode**  
   - When scanning many `k`, use a Sieve of Eratosthenes to precompute all primes up to `8N - 1`.  
   - Avoids repeated Miller–Rabin calls.

---

## Usage Examples
Just look at the code in `main()`, but here are the same examples.

Depending on your value of `k`, you will need to up the `m_limit` parameter. The higher the `m_limit` value, the longer it is going to take but we need to be able to up it for certain massive values of `k`.


### Generate all valid values for `k = 1..N`

```python
from solver import generate_m_map

# Generate mapping k -> m for k up to 200
m_map = generate_m_map(200, m_limit=500)

# Print nicely
for k, m in sorted(m_map.items()):
    print(f"k = {k}, m = {m}")
```

### Example Output:
```
k = 3, m = 2
k = 4, m = 2
k = 6, m = 6
k = 9, m = 6
k = 10, m = 8
...
```

### Find `m` for a single `k`
```
result = find_m_for_k(4, m_limit=1000)
```

### Run sanity checks
```
m_map = generate_m_map(50, m_limit=500)
sanity_check_map(m_map)  # prints only if something is wrong
```