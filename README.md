# Scalar Multiplication on Elliptic Curves over Finite Fields in Rust

## Overview

This Rust project implements scalar multiplication on elliptic curves defined over finite fields $\mathbb{Z}_p$. The primary focus is to compute $nP$ for a given point $P$ using elliptic curve arithmetic. This is a crucial operation in elliptic curve cryptography (ECC), which is widely used in secure communication, digital signatures, and key exchange protocols.

### Key Features
- **Elliptic Curve Point Addition**: Handles the addition of two points on the elliptic curve.
- **Point Doubling**: Computes the result of doubling a point.
- **Scalar Multiplication**: Efficiently computes $nP$ using a double-and-add algorithm.
- **Curve Validation**: Checks whether a point lies on the given elliptic curve.

## Mathematical Background

The elliptic curve used is defined by the equation:

$y^2 \equiv x^3 + ax + b \pmod{p}$
where:
- $a$ and $b$ are coefficients defining the curve.
- $p$ is a prime number representing the finite field $\mathbb{Z}_p$.
- A point $P = (x, y)$ satisfies this equation if it lies on the curve.

### Point Addition and Doubling

Given two points $P = (x_1, y_1)$ and $Q = (x_2, y_2)$ on the curve:
1. **Point Addition** ($P \neq Q$):

   $\lambda = \frac{y_2 - y_1}{x_2 - x_1} \pmod{p}$
   
   $x_3 = \lambda^2 - x_1 - x_2 \pmod{p},$

   $\quad y_3 = \lambda (x_1 - x_3) - y_1 \pmod{p}$
   

2. **Point Doubling** ($P = Q $):
   
   $\lambda = \frac{3x_1^2 + a}{2y_1} \pmod{p}$
   
   $x_3 = \lambda^2 - 2x_1 \pmod{p},$

   $y_3 = \lambda (x_1 - x_3) - y_1 \pmod{p}$

3. **Scalar Multiplication** ($nP$):
   Uses the double-and-add method to compute $nP$ efficiently.

## Implementation

### Project Structure
- **`Point` Struct**: Represents a point on the elliptic curve, including support for a point at infinity.
- **`mod_inv` Function**: Computes the modular inverse using the Extended Euclidean Algorithm.
- **`elliptic_add` Function**: Adds two points on the curve.
- **`scalar_mult` Function**: Performs scalar multiplication using a double-and-add approach.
- **`is_on_curve` Function**: Checks whether a point is on the elliptic curve.

### Code Explanation

```rust
#[derive(Clone, Debug)]
struct Point {
    x: i64,
    y: i64,
    infinity: bool,
}

impl Point {
    fn is_at_infinity(&self) -> bool {
        self.infinity
    }

    fn at_infinity() -> Self {
        Point { x: 0, y: 0, infinity: true }
    }
}

// Function to compute the modular inverse
fn mod_inv(a: i64, m: i64) -> i64 {
    let (mut t, mut new_t, mut r, mut new_r) = (0, 1, m, a % m);
    while new_r != 0 {
        let quotient = r / new_r;
        t = t - quotient * new_t;
        std::mem::swap(&mut t, &mut new_t);
        r = r - quotient * new_r;
        std::mem::swap(&mut r, &mut new_r);
    }
    if r > 1 {
        panic!("{} has no modular inverse modulo {}", a, m);
    }
    if t < 0 { t += m; }
    t
}

// Elliptic curve point addition
fn elliptic_add(p: &Point, q: &Point, a: i64, m: i64) -> Point {
    if p.is_at_infinity() {
        return q.clone();
    }
    if q.is_at_infinity() {
        return p.clone();
    }
    if p.x == q.x && (p.y != q.y || p.y == 0) {
        return Point::at_infinity();
    }

    let (x1, y1) = (p.x, p.y);
    let (x2, y2) = (q.x, q.y);

    let lambda = if x1 == x2 && y1 == y2 {
        let num = (3 * x1 * x1 + a).rem_euclid(m);
        let denom = mod_inv(2 * y1, m);
        (num * denom).rem_euclid(m)
    } else {
        let num = (y2 - y1).rem_euclid(m);
        let denom = mod_inv((x2 - x1).rem_euclid(m), m);
        (num * denom).rem_euclid(m)
    };

    let x3 = (lambda * lambda - x1 - x2).rem_euclid(m);
    let y3 = (lambda * (x1 - x3) - y1).rem_euclid(m);

    Point {
        x: if x3 < 0 { x3 + m } else { x3 },
        y: if y3 < 0 { y3 + m } else { y3 },
        infinity: false,
    }
}

// Scalar multiplication
fn scalar_mult(n: i64, p: &Point, a: i64, m: i64) -> Point {
    let mut result = Point::at_infinity();
    let mut addend = p.clone();
    let mut n = n;

    while n > 0 {
        if n & 1 == 1 {
            result = elliptic_add(&result, &addend, a, m);
        }
        addend = elliptic_add(&addend, &addend, a, m);
        n >>= 1;
    }

    result
}

fn main() {
    let p = Point { x: 205, y: 130, infinity: false };
    let n = 2;
    let a = 4;
    let b = 4;
    let m = 313;

    // Compute nP
    let np = scalar_mult(n, &p, a, m);
    if np.is_at_infinity() {
        println!("{}P is the point at infinity", n);
    } else {
        println!("{}P = ({}, {})", n, np.x, np.y);
        println!("Is the point on the curve? {}", is_on_curve(&np, a, b, m));
    }
}
```
## Requirements
- To get started, ensure you have [Rust](https://www.rust-lang.org/tools/install) installed on your machine. You can then clone the repository and build the project.

## Installation

1. **Clone the repository**:
    ```bash
    git clone https://github.com/cypriansakwa/Elliptic_Curves_Scalar_Multiplication_over_Finite_Fields_in_Rust.git
    cd Elliptic_Curves_Scalar_Multiplication_over_Finite_Fields_in_Rust
    ```

2. **Build the project**:
    ```bash
    cargo build
    ```

3. **Run the program**:
    ```bash
    cargo run
    ```

## Usage

To execute the program, simply run:

```bash
cargo run
```
