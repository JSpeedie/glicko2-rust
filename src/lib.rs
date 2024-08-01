#![warn(missing_docs)]

//! # Introduction
//!
//! To fully understand the functions contained in this library it will be
//! necessary to understand the names Glicko-2 uses for various variables.
//! Glicko-2 uses 3 variables to represent a given player in the system:
//!
//! 1. The player's rating (*r*)
//! 2. The player's rating deviation (RD)
//! 3. The player's rating volatility (σ).
//!
//! Glicko-2 also uses 1 system constant (τ) in some of the calculations of
//! the algorithm. What is important to note here is that two of these variables
//! are stored in one scale or range, and used for calculations in another
//! scale. The player's rating is stored in the "original scale" (i.e. the
//! original Glicko scale) and is converted to the "Glicko-2 scale" for
//! calculations. Likewise, the player's rating deviation is stored in the
//! "original scale" and converted to the "Glicko-2 scale" for calculations.
//! You can refer to the following table to understand which scale a given
//! variable is on based on the name that is being used to refer to it:
//!
//! | Original Scale | Glicko-2 Scale |
//! | -------------- | -------------- |
//! | *r*            | μ              |
//! | RD             | Φ              |
//!
//! It should also be noted that the volatility (σ) does not change between
//! scales and neither does the system constant (τ).


use std::f64::consts::E;
use std::f64::consts::PI;


struct G2Context {
    tau: f64, // The system constant which "constrains the change in volatility over time"
}

pub struct G2Player {
    rating: f64,
    rd: f64, // The rating deviation of the player
    /// Reflects "the degree of expected fluctuation in a player's rating"
    pub vol: f64,
}

impl G2Player {
    /// Returns the value of a given `G2Player`'s rating (μ) on the Glicko-2
    /// scale. This function paired with the `get_rd()` function constitute
    /// Step 2 of the Glicko-2 algorithm.
    ///
    /// Parameters:
    /// * 'self' a reference to the `G2Player` struct whose Glicko-2 scale
    ///     rating will be returned.
    ///
    /// Return:
    /// * An f64 representing the player's rating converted to the Glicko-2
    ///     scale.
    pub fn get_rating(&self) -> f64 {
        (self.rating * 173.7178) + 1500.0
    }

    /// Sets the `G2Player`'s rating (*r*) on the original scale based on the
    /// given Glicko-2 scale value `new_rating`. This function paired with the
    /// `set_rd()` function consitute Step 8 of the Glicko-2 Algorithm.
    ///
    /// Parameters:
    /// * `self` a reference to the `G2Player` struct whose rating will be
    ///     changed.
    /// * `new_rating` the Glicko-2 scale value to set `self`'s rating to.
    ///
    /// Return:
    /// * Nothing - the function modifies the `G2Player` struct.
    pub fn set_rating(&mut self, new_rating: f64) {
        self.rating = (new_rating - 1500.0) / 173.7178;
    }

    /// Returns the value of a given `G2Player`'s rating deviation (Φ) on the
    /// Glicko-2 scale. This function paired with the `get_rating()` function
    /// constitute Step 2 of the Glicko-2 algorithm.
    ///
    /// Parameters:
    /// * 'self' a reference to the `G2Player` struct whose Glicko-2 scale
    ///     rating deviation will be returned.
    ///
    /// Return:
    /// * An f64 representing the player's rating deviation converted to the
    ///     Glicko-2 scale.
    pub fn get_rd(&self) -> f64 {
        self.rd * 173.7178
    }

    /// Sets the `G2Player`'s rating deviation (RD) on the original scale based on
    /// the given Glicko-2 scale value `new_rd`. This function paired with the
    /// `set_rating()` function consitute Step 8 of the Glicko-2 Algorithm.
    ///
    /// Parameters:
    /// * `self` a reference to the `G2Player` struct whose rating deviation
    ///     will be changed.
    /// * 'new_rd' the Glicko-2 scale value to set `self`'s rating deviation to.
    ///
    /// Return:
    /// * Nothing - the function modifies the `G2Player` struct.
    pub fn set_rd(&mut self, new_rd: f64) {
        self.rd = new_rd / 173.7178;
    }
}

/// The Glicko-2 *g* function used in Steps 3, 4 and 7.
///
/// Parameters:
/// * `phi` the Glicko-2 scale rating deviation of the player we need to call
///     this function for.
///
/// Return:
/// * The exact purpose or meaning of the value returned by this function is
///     not made clear by the Glicko-2 PDF.
pub fn _g(phi: f64) -> f64 {
	1.0 / (1.0 + (3.0 * phi.powi(2) / PI.powi(2))).sqrt()
}


/// The Glicko-2 *E* function used in Steps 3, 4, and 7.
///
/// Parameters:
/// * `mu` the Glicko-2 scale rating of the player we need to call
///     this function for.
/// * `mu_j` the Glicko-2 scale rating of an opponent the player-of-interest
///     played in the last rating period.
/// * `phi_j` the Glicko-2 scale rating deviation of an opponent the
///     player-of-interest played in the last rating period.
///
/// Return:
/// * The exact purpose or meaning of the value returned by this function is
///     not made clear by the Glicko-2 PDF.
pub fn _E(mu: f64, mu_j: f64, phi_j: f64) -> f64 {
	1.0 / (1.0 + E.powf(-1.0 * _g(phi_j) * (mu - mu_j)))
}


/// Calculate the v (*v*) quantity representing the estimated variance in
/// a player's rating based only on game outcomes. This function constitutes
/// Step 3 of the Glicko-2 algorithm.
///
/// Parameters:
/// * `mu` the Glicko-2 scale rating of the player we need to call
///     this function for.
/// * `opp_ratings` an array of f64s containing (in order) all the ratings of
///     the opponents `player` played.
/// * `opp_rds` an array of f64s containing (in order) all the RDs of the
///     opponents `player` played.
///
/// Return:
/// * The v (*v*) quantity representing the estimated variance in
///     a player's rating based only on game outcomes.
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
/// Parameters:
/// * `mu` the Glicko-2 scale rating of the player we need to call
///     this function for.
/// * `opp_ratings` an array of f64s containing (in order) all the ratings of
///     the opponents `player` played.
/// * `opp_rds` an array of f64s containing (in order) all the RDs of the
///     opponents `player` played.
/// * `scores` an array of f64s containing (in order) the scores
///     of the games `player` played. 1 for `player` winning, 0.5 for a draw and
///     0 for a `player` loss
/// * `v` the quantity *v* calculated in Step 3 from this series of outcomes
///     representing the estimated variance of the player.
///
/// Return:
/// * The estimated improvement in rating for a player who has the played
///     the given series of games.
pub fn delta(mu: f64, opp_ratings: &[f64], opp_rds: &[f64], scores: &[f64],
             v: f64) -> f64 {

    let mut sum: f64 = 0.0;

    for j in 0..opp_ratings.len() {
        sum += _g(opp_rds[j]) * (scores[j] - _E(mu, opp_ratings[j], opp_rds[j]));
    }

    v * sum
}


/// Calculates the new volatility (σ') of the player with the rating deviation
/// `phi` (Φ), and the volatility (σ) `sigma`, in a system with a system
/// constant (τ) `tau`, after a series of outcomes with a delta (∆) `delta` and
/// a *v* `v`. This function constitutes Step 5 of the Glicko-2 algorithm.
///
/// * `phi` the Glicko-2 scale value representing the player's rating deviation.
/// * `sigma` the value representing the player's volatility.
/// * `tau` the system constant.
/// * `v` the quantity *v* calculated in Step 3 from a series of outcomes
///     representing the estimated variance of the player.
/// * `delta` the quantity ∆ calculated in Step 4 from a series of outcomes
///     representing the estimated improvement in rating of the player.
///
/// Return:
/// * An f64 representing the player's new volatility (σ').
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
    let mut fA = f(A);
    let mut fB = f(B);
    let mut fC;

    // Step 5, item 4
    while (B - A).abs() > epsilon {
        // item 4a
        let C = A + (((A - B) * fA) / (fB - fA));
        fC = f(C);
        // item 4b
        if (fC * fB) <= 0.0 {
            A = B;
            fA = fB;
        } else {
            fA = fA / 2.0;
        }
        // item 4c
        B = C;
        fB = fC;
    }
    // Step 5, item 5
    E.powf(A / 2.0)
}


/// Returns a provisional Glicko-2 scale rating deviation (Φ*) representing the
/// new pre-rating-period value. This function constitutes Step 6 of the
/// Glicko-2 algorithm.
///
/// This function can be used to increase the RD of players who did not compete
/// in the rating period by providing their current rating deviation and
/// current volatility (Φ, σ) instead of their current rating deviation
/// and their new volatility (Φ, σ').
///
/// Parameters:
/// * `phi` the Glicko-2 scale rating deviation of the player we are performing
///     this step of the algorithm for.
/// * `new_sigma` the new volatility of the player we are performing this step
///     of the algorithm for, calculated in Step 5.
///
/// Return:
/// * An f64 representing the player's provisional rating deviation (Φ*)
///     calculated in Step 6.
pub fn update_rd(phi: f64, new_sigma: f64) -> f64 {
    (phi.powi(2) + new_sigma.powi(2)).sqrt()
}


/// Returns an updated rating (μ') and an updated rating deviation (Φ') for a
/// player that undergoes a series of outcomes with the given scores against
/// opponents with the given ratings and RDs. This function constitutes Step 7
/// of the Glicko-2 algorithm.
///
/// Parameters:
/// * `mu` the Glicko-2 scale rating of the player we need to call this
///     function for.
/// * `prov_phi` the Glicko-2 scale provisional rating deviation of the player
///     we are performing this step of the algorithm for. This value is
///     calculated by the `update_rd()` function aka. Step 6.
/// * 'opp_ratings' an array of f64s containing (in order) all the ratings of
///     the opponents `player` played.
/// * 'opp_rds' an array of f64s containing (in order) all the RDs of the
///     opponents `player` played.
/// * 'scores' an array of f32s containing (in order) the scores
///     of the games `player` played. 1 for `player` winning, 0.5 for a draw and
///     0 for a `player` loss
///
/// Return:
/// * A tuple containing the player's new rating and RD on the Glicko-2 scale
///     (Φ', μ')
pub fn update_player(mu: f64, prov_phi: f64, opp_ratings: &[f64],
                     opp_rds: &[f64], scores: &[f64], v: f64) -> (f64, f64) {

    let new_phi = 1.0 / ((1.0 / prov_phi.powi(2)) + (1.0 / v)).sqrt();
    let mut new_mu: f64 = 0.0;

    for j in (0..opp_ratings.len()) {
        new_mu += _g(opp_rds[j]) * (scores[j] - _E(mu, opp_ratings[j], opp_rds[j]));
    }
    new_mu = (mu + new_phi.powi(2)) * new_mu;

    (new_phi, new_mu)
}


// #[cfg(test)]
// mod tests {
//     use super::*;
// 
//     #[test]
//     fn it_works() {
//         let result = add(2, 2);
//         assert_eq!(result, 4);
//     }
// }
