#[derive(Clone, Debug)]
struct Point {
    x: i64,
    y: i64,
    infinity: bool,  // Flag to indicate if the point is the point at infinity
}

impl Point {
    fn is_at_infinity(&self) -> bool {
        self.infinity
    }

    fn at_infinity() -> Self {
        Point { x: 0, y: 0, infinity: true }
    }
}

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
        // Point doubling
        let num = (3 * x1 * x1 + a).rem_euclid(m);
        let denom = mod_inv(2 * y1, m);
        (num * denom).rem_euclid(m)
    } else {
        // Point addition
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

fn is_on_curve(p: &Point, a: i64, b: i64, m: i64) -> bool {
    if p.is_at_infinity() {
        return true;
    }
    let left_side = (p.y * p.y).rem_euclid(m);
    let right_side = ((p.x * p.x * p.x) + a * p.x + b).rem_euclid(m);
    left_side == right_side
}

fn main() {
    let p = Point { x: 205, y: 130, infinity: false };
    let n = 2;
    let a = 4;
    let b = 4;
    let m = 313;

    // Compute and print the final result
    let np = scalar_mult(n, &p, a, m);
    if np.is_at_infinity() {
        println!("{}P is the point at infinity", n);
    } else {
        println!("{}P = ({}, {})", n, np.x, np.y);
        println!("Is the point on the curve? {}", is_on_curve(&np, a, b, m));
    }
}












