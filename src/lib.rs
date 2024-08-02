#![warn(missing_docs)]

//! # Preface
//!
//! To fully understand the functions contained in this library it will be
//! necessary to understand the names Glicko-2 uses for various variables.
//! Glicko-2 uses 3 variables to represent a given player in the system:
//!
//! 1. The player's rating (*r*)
//! 2. The player's rating deviation (RD)
//! 3. The player's rating volatility (σ).
//!
//! Glicko-2 also uses one system constant (τ) in some of the calculations of
//! the algorithm. What is important to note here is that two of these variables
//! are stored in one scale or range, and then converted to another scale for
//! calculations. The player's rating is stored in the "original scale" (i.e.
//! the original Glicko scale) and is converted to the "Glicko-2 scale" for
//! calculations. Likewise, the player's rating deviation is stored in the
//! "original scale" and converted to the "Glicko-2 scale" for calculations.
//! You can refer to the following table to understand which scale a given
//! variable is on based on the name that is being used to refer to it:
//!
//! | Variable         | Original Scale | Glicko-2 Scale |
//! | ---------------- | -------------- | -------------- |
//! | Rating           | *r*            | μ              |
//! | Rating Deviation | RD             | Φ              |
//!
//! It should also be noted that the volatility (σ) does not change between
//! scales and neither does the system constant (τ).


use std::f64::consts::{E, PI};


/// Returns the given original scale rating (*r*) converted to the Glicko-2 scale (μ). This
/// function paired with the `rd_to_g2_scale()` function constitute Step 2 of the Glicko-2 algorithm.
///
/// #### Parameters:
/// * `rating` the original scale rating (*r*) to be converted.
/// #### Return:
/// * An `f64` representing the rating converted to the Glicko-2 scale (μ).
pub fn rating_to_g2_scale(rating: f64) -> f64 {
    (rating - 1500.0) / 173.7178
}


/// Returns the given original scale rating deviation (RD) converted to the Glicko-2 scale (Φ).
/// This function paired with the `rating_to_g2_scale()` function constitute Step 2 of the Glicko-2
/// algorithm.
///
/// #### Parameters:
/// * `rd` the original scale rating deviation (RD) to be converted.
/// #### Return:
/// * An `f64` representing the rating deviation converted to the Glicko-2 scale (Φ).
pub fn rd_to_g2_scale(rd: f64) -> f64 {
    rd / 173.7178
}


/// The Glicko-2 *g* function used in Steps 3, 4, 7, and by the *E* function.
///
/// #### Parameters:
/// * `phi` the Glicko-2 scale rating deviation (Φ) of the player we need to
///     call this function for.
/// #### Return:
/// * An `f64` produced by the *g* function of the Glicko-2 algorithm. The
///     exact purpose or meaning of the value returned by this function is not
///     made clear by the Glicko-2 PDF.
pub fn _g(phi: f64) -> f64 {
	1.0 / (1.0 + (3.0 * phi.powi(2) / PI.powi(2))).sqrt()
}


/// The Glicko-2 *E* function used in Steps 3, 4, and 7.
///
/// #### Parameters:
/// * `mu` the Glicko-2 scale rating (μ) of the player we need to call this
///     function for.
/// * `mu_j` the Glicko-2 scale rating (μ<sub>j</sub>) of the player-of-interest's
///     j<sup>th</sup> opponent played in the last rating period.
/// * `phi_j` the Glicko-2 scale rating deviation (Φ<sub>j</sub>) of the
///     player-of-interest's j<sup>th</sup> opponent played in the last rating
///     period.
/// #### Return:
/// * An `f64` produced by the *E* function of the Glicko-2 algorithm. The
///     exact purpose or meaning of the value returned by this function is not
///     made clear by the Glicko-2 PDF.
#[allow(non_snake_case)]
pub fn _E(mu: f64, mu_j: f64, phi_j: f64) -> f64 {
	1.0 / (1.0 + E.powf(-1.0 * _g(phi_j) * (mu - mu_j)))
}


/// Calculate the v (*v*) quantity representing the estimated variance in
/// a player's rating based only on game outcomes. This function constitutes
/// Step 3 of the Glicko-2 algorithm.
///
/// #### Parameters:
/// * `mu` the Glicko-2 scale rating of the player we need to call
///     this function for.
/// * `opp_ratings` an array (μ<sub>1</sub>, ... , μ<sub>m</sub>) of `f64`s
///     containing (in order) all the ratings of the opponents the
///     player-of-interest played.
/// * `opp_rds` an array (Φ<sub>1</sub>, ... , Φ<sub>m</sub>) of `f64`s
///     containing (in order) all the rating deviations of the opponents
///     the player-of-interest played.
/// #### Return:
/// * An `f64` representing the estimated variance in a player's rating (*v*)
///     for a player who has played the given series of games.
pub fn _v(mu: f64, opp_ratings: &[f64], opp_rds: &[f64]) -> f64 {
    let mut sum: f64 = 0.0;

    for j in 0..opp_ratings.len() {
        sum +=
            _g(opp_rds[j]).powi(2)
            * _E(mu, opp_ratings[j], opp_rds[j])
            * (1.0 - _E(mu, opp_ratings[j], opp_rds[j]));
    }

    sum.powi(-1)
}


/// Calculate the delta (∆) quantity representing the estimated improvement
/// in rating based only on game outcomes. This function constitutes Step 4 of
/// the Glicko-2 algorithm.
///
/// #### Parameters:
/// * `mu` the Glicko-2 scale rating (μ) of the player we need to call this
///     function for.
/// * `opp_ratings` an array (μ<sub>1</sub>, ... , μ<sub>m</sub>) of `f64`s
///     containing (in order) all the ratings of the opponents the
///     player-of-interest played.
/// * `opp_rds` an array (Φ<sub>1</sub>, ... , Φ<sub>m</sub>) of `f64`s
///     containing (in order) all the rating deviations of the opponents
///     the player-of-interest played.
/// * `scores` an array (s<sub>1</sub>, ... , s<sub>m</sub>) of `f64`s
///     containing (in order) the scores of the games the player-of-interest
///     played. A value of 1 indicates that the player-of-interest won, 0.5
///     indicates a draw, and 0 indicates that the player-of-interest lost.
/// * `v` the quantity *v* calculated in Step 3 from this series of outcomes
///     representing the estimated variance of the player.
/// #### Return:
/// * An `f64` representing the estimated improvement in rating (∆) for a player
///     who has the played the given series of games.
pub fn delta(mu: f64, opp_ratings: &[f64], opp_rds: &[f64], scores: &[f64],
             v: f64) -> f64 {

    let mut sum: f64 = 0.0;

    for j in 0..opp_ratings.len() {
        sum += _g(opp_rds[j]) * (scores[j] - _E(mu, opp_ratings[j], opp_rds[j]));
    }

    v * sum
}


/// Calculates the new volatility (σ') of the player with the rating deviation
/// (Φ) `phi`, and the volatility (σ) `sigma`, in a system with a system
/// constant (τ) `tau`, after a series of outcomes with a delta (∆) `delta` and
/// a *v* `v`. This function constitutes Step 5 of the Glicko-2 algorithm.
///
/// #### Parameters:
/// * `phi` the Glicko-2 scale value representing the player's rating deviation
///     (Φ).
/// * `sigma` the value representing the player's volatility (σ).
/// * `tau` the system constant (τ).
/// * `v` the quantity *v* calculated in Step 3 from a series of outcomes
///     representing the estimated variance of the player.
/// * `delta` the quantity ∆ calculated in Step 4 from a series of outcomes
///     representing the estimated improvement in rating of the player.
/// #### Return:
/// * An `f64` representing the player's new volatility (σ').
#[allow(non_snake_case)]
pub fn new_vol(phi: f64, sigma: f64, tau: f64, v: f64, delta: f64) -> f64 {
    // Step 5, item 1
    let a: f64 = (sigma.powi(2)).ln();

    let f = |x: f64| -> f64 {
        let numerator1 = E.powf(x) * (delta.powi(2) - phi.powi(2) - v - E.powf(x));
        let denominator1 = 2.0 * (phi.powi(2) + v + E.powf(x)).powi(2);

        (numerator1 / denominator1) - ((x - a) / tau.powi(2))
    };
    let epsilon = 0.000001;

    // Step 5, item 2
    let mut A = a;
    let mut B;
    if delta.powi(2) > phi.powi(2) + v {
        B = (delta.powi(2) - phi.powi(2) - v).ln();
    } else {
        let mut k: f64 = 1.0;
        while f(a - (k * tau)) < 0.0 {
            k = k + 1.0;
        }
        B = a - (k * tau);
    }

    // Step 5, item 3
    let mut f_A = f(A);
    let mut f_B = f(B);
    let mut f_C;

    // Step 5, item 4
    while (B - A).abs() > epsilon {
        // item 4a
        let C = A + (((A - B) * f_A) / (f_B - f_A));
        f_C = f(C);
        // item 4b
        if (f_C * f_B) <= 0.0 {
            A = B;
            f_A = f_B;
        } else {
            f_A /= 2.0;
        }
        // item 4c
        B = C;
        f_B = f_C;
    }
    // Step 5, item 5
    E.powf(A / 2.0)
}


/// Returns a provisional Glicko-2 scale rating deviation (Φ*) representing the
/// new pre-rating-period value. This function constitutes Step 6 of the
/// Glicko-2 algorithm.
///
/// This function can be used to increase the rating deviation of players who
/// did not compete in the rating period by calling this function with their
/// current rating deviation and current volatility (Φ, σ) instead of performing
/// Steps 3-7 of the Glicko-2 algorithm. For reference, for players who did
/// compete, this function is called in Step 6 with their current rating
/// deviation and their new volatility (Φ, σ').
///
/// #### Parameters:
/// * `phi` the Glicko-2 scale rating deviation (Φ) of the player we are
///     performing this step of the algorithm for.
/// * `new_sigma` the new volatility (σ') of the player we are performing this
///     step of the algorithm for, calculated in Step 5.
/// #### Return:
/// * An `f64` representing the player's provisional rating deviation (Φ*)
///     calculated in Step 6.
pub fn update_rd(phi: f64, new_sigma: f64) -> f64 {
    (phi.powi(2) + new_sigma.powi(2)).sqrt()
}


/// Returns an updated rating (μ') and an updated rating deviation (Φ') for a
/// player that undergoes a series of outcomes with scores
/// (s<sub>1</sub>, ... , s<sub>m</sub>) against opponents with the ratings
/// (μ<sub>1</sub>, ... , μ<sub>m</sub>) and rating deviations
/// (Φ<sub>1</sub>, ... , Φ<sub>m</sub>). This function constitutes Step 7 of
/// the Glicko-2 algorithm.
///
/// #### Parameters:
/// * `mu` the Glicko-2 scale rating (μ) of the player we need to call this
///     function for.
/// * `prov_phi` the Glicko-2 scale provisional rating deviation (Φ*) of the
///     player we are performing this step of the algorithm for. This value is
///     calculated by the `update_rd()` function aka. Step 6.
/// * `opp_ratings` an array (μ<sub>1</sub>, ... , μ<sub>m</sub>) of `f64`s
///     containing (in order) all the ratings of the opponents the
///     player-of-interest played.
/// * `opp_rds` an array (Φ<sub>1</sub>, ... , Φ<sub>m</sub>) of `f64`s
///     containing (in order) all the rating deviations of the opponents
///     the player-of-interest played.
/// * `scores` an array (s<sub>1</sub>, ... , s<sub>m</sub>) of `f64`s
///     containing (in order) the scores of the games the player-of-interest
///     played. A value of 1 indicates that the player-of-interest won, 0.5
///     indicates a draw, and 0 indicates that the player-of-interest lost.
/// #### Return:
/// * An `(f64, f64)` tuple containing the player-of-interest's new rating and
///     rating deviation on the Glicko-2 scale (Φ', μ').
pub fn update_player(mu: f64, prov_phi: f64, opp_ratings: &[f64],
                     opp_rds: &[f64], scores: &[f64], v: f64) -> (f64, f64) {

    let new_phi = 1.0 / ((1.0 / prov_phi.powi(2)) + (1.0 / v)).sqrt();
    let mut new_mu: f64 = 0.0;

    for j in 0..opp_ratings.len() {
        new_mu += _g(opp_rds[j]) * (scores[j] - _E(mu, opp_ratings[j], opp_rds[j]));
    }
    new_mu = (mu + new_phi.powi(2)) * new_mu;

    (new_phi, new_mu)
}


/// Returns the given Glicko-2 scale rating (μ) converted to the original scale (*r*). This
/// function paired with the `rd_to_original_scale()` function constitute Step 8 of the Glicko-2
/// algorithm.
///
/// #### Parameters:
/// * `rating` the Glicko-2 scale rating (μ) to be converted.
/// #### Return:
/// * An `f64` representing the rating converted to the original scale (*r*).
pub fn rating_to_original_scale(rating: f64) -> f64 {
    (rating * 173.7178) + 1500.0
}


/// Returns the given Glicko-2 scale rating deviation (Φ) converted to the original scale (RD).
/// This function paired with the `rating_to_original_scale()` function constitute Step 8 of the
/// Glicko-2 algorithm.
///
/// #### Parameters:
/// * `rd` the Glicko-2 scale rating deviation (Φ) to be converted.
/// #### Return:
/// * An `f64` representing the rating deviation converted to the original scale (RD).
pub fn rd_to_original_scale(rd: f64) -> f64 {
    rd * 173.7178
}


#[cfg(test)]
mod tests {
    use super::*; // Use all items from parent module
    use assert_float_eq::*; // Get float comparison functions

    #[test]
    fn scale_rating_to_g2_scale() {
        let ratings: [f64; 4] = [1500.0, 1400.0, 1550.0, 1700.0];

        assert_f64_near!(rating_to_g2_scale(ratings[0]), 0.0);
        assert_f64_near!(rating_to_g2_scale(ratings[1]), -0.5756462492617337);
        assert_f64_near!(rating_to_g2_scale(ratings[2]), 0.28782312463086684);
        assert_f64_near!(rating_to_g2_scale(ratings[3]), 1.1512924985234674);
    }

    #[test]
    fn scale_rd_to_g2_scale() {
        let rds: [f64; 4] = [200.0, 30.0, 100.0, 300.0];

        assert_f64_near!(rd_to_g2_scale(rds[0]), 1.1512924985234674);
        assert_f64_near!(rd_to_g2_scale(rds[1]), 0.1726938747785201);
        assert_f64_near!(rd_to_g2_scale(rds[2]), 0.5756462492617337);
        assert_f64_near!(rd_to_g2_scale(rds[3]), 1.726938747785201);
    }

    #[test]
    fn g_function() {
        let phis: [f64; 4] = [1.1512924985234674, 0.1726938747785201,
                              0.5756462492617337, 1.726938747785201];

        assert_f64_near!(_g(phis[0]), 0.8442815045053368);
        assert_f64_near!(_g(phis[1]), 0.9954980064506083);
        assert_f64_near!(_g(phis[2]), 0.9531489778689763);
        assert_f64_near!(_g(phis[3]), 0.7242354780877525);
    }

    #[test]
    fn e_function() {
        let mu: f64 = 0.0;
        let mus: [f64; 3] = [-0.5756462492617337, 0.28782312463086684, 1.1512924985234674];
        let phis: [f64; 3] = [0.1726938747785201, 0.5756462492617337, 1.726938747785201];

        assert_f64_near!(_E(mu, mus[0], phis[0]), 0.6394677305521533);
        assert_f64_near!(_E(mu, mus[1], phis[1]), 0.4318423561076679);
        assert_f64_near!(_E(mu, mus[2], phis[2]), 0.3028407290952193);
    }
}
