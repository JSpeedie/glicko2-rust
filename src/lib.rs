/// Introduction:
///
/// To fully understand the functions contained in this library it will be
/// necessary to understand the names Glicko-2 uses for various variables.
/// Glicko-2 uses 3 variables to represent a given player in the system:
///
/// 1. The player's rating (*r*)
/// 2. The player's rating deviation (RD)
/// 3. The player's rating volatility (σ).
///
/// Glicko-2 also uses 1 system constant (τ) in some of the calculations of
/// the algorithm. What is important to note here is that two of these variables
/// are stored in one scale or range, and used for calculations in another
/// scale. The player's rating is stored in the "original scale" (i.e. the
/// original Glicko scale) and is converted to the "Glicko-2 scale" for
/// calculations. Likewise, the player's rating deviation is stored in the
/// "original scale" and converted to the "Glicko-2 scale" for calculations.
/// You can refer to the following table to understand which scale a given
/// variable is on based on the name that is being used to refer to it:
///
/// | Original Scale         | Glicko-2 Scale |
/// | ---------------------- | -------------- |
/// | *r* or rating          | μ              |
/// | RD or rating deviation | Φ              |
///
/// It should also be noted that the volatility (σ) does not change between
/// scales and neither does the system constant (τ).


use std::f64::consts::E;

struct G2Context {
    tau: f64, // The system constant which "constrains the change in volatility over time"
}

struct G2Player {
    rating: f64,
    rd: f64, // The rating deviation of the player
    vol: f64, // Reflects "the degree of expected fluctuation in a player's rating"
}


/// Returns the value of a given `G2Player`'s rating (μ) on the Glicko-2 scale.
/// This function paired with the `getRD()` function constitute Step 2 of the
/// Glicko-2 algorithm.
///
/// Parameters:
/// * 'player' a reference to the `G2Player` struct whose Glicko-2 scale rating
///     will be returned.
/// Return:
/// * An f64 representing the player's rating converted to the Glicko-2 scale.
fn getRating(player: &G2Player) -> f64 {
    (player.rating * 173.7178) + 1500.0
}


/// Sets the `G2Player`'s rating (*r*) on the original scale based on the given
/// Glicko-2 scale value `new_rating`. This function paired with the `setRD()`
/// function consitute Step 8 of the Glicko-2 Algorithm.
///
/// Parameters:
/// * `player` the `G2Player` struct whose rating will be changed.
/// * `new_rating` the Glicko-2 scale value to set `player`'s rating to.
/// Return:
/// * Nothing - the function modifies the `G2Player` struct.
fn setRating(player: &mut G2Player, new_rating: f64) {
	player.rating = (new_rating - 1500.0) / 173.7178;
}


/// Returns the value of a given `G2Player`'s rating deviation (Φ) on the
/// Glicko-2 scale. This function paired with the `getRating()` function
/// constitute Step 2 of the Glicko-2 algorithm.
///
/// Parameters:
/// * 'player' a reference to the `G2Player` struct whose Glicko-2 scale rating
///     deviation will be returned.
/// Return:
/// * An f64 representing the player's rating deviation converted to the
///     Glicko-2 scale.
fn getRd(player: &G2Player) -> f64 {
	player.rd * 173.7178
}


/// Sets the `G2Player`'s rating deviation (RD) on the original scale based on
/// the given Glicko-2 scale value `new_rd`. This function paired with the
/// `setRating()` function consitute Step 8 of the Glicko-2 Algorithm.
///
/// Parameters:
/// * `player` the `G2Player` struct whose rating deviation will be changed.
/// * 'new_rd' the Glicko-2 scale value to set `player`'s rating deviation to.
/// Return:
/// * Nothing - the function modifies the `G2Player` struct.
fn setRd(player: &mut G2Player, new_rd: f64) {
	player.rd = new_rd / 173.7178;
}


/// Returns a provisional Glicko-2 scale rating deviation (Φ*) representing the
/// new pre-rating-period value. This function constitutes Step 6 of the
/// Glicko-2 algorithm.
///
/// Parameters:
/// * `phi` the Glicko-2 scale rating deviation of the player we are performing
///     this step of the algorithm for.
/// * `sigma` the volatility of the player we are performing this step of the
///     algorithm for.
/// Return:
/// * An f64 representing the player's provisional rating deviation calculated
///     in Step 6.
fn updateRD(phi: f64, sigma: f64) {
    (phi.powi(2) + sigma.powi(2)).sqrt();
}


/// Updates a `G2Player`'s rating after the player undergoes a series of
/// outcomes with the given scores against opponents with the given ratings
/// and RDs.
///
/// Parameters:
/// * `player` a `G2Player` struct representing the player whose rating is to
///     be updated.
/// * 'opp_ratings' an array of f64s containing (in order) all the ratings of
///     the opponents `player` played.
/// * 'opp_rds' an array of f64s containing (in order) all the RDs of the
///     opponents `player` played.
/// * 'scores' an array of f32s containing (in order) the scores
///     of the games `player` played. 1 for `player` winning, 0.5 for a draw and
///     0 for a `player` loss
/// Return:
/// * Nothing - the function modifies the `G2Player` struct.
fn update_player(player: &mut G2Player, opp_ratings: &f64, opp_rds: &f64,
                 scores: &f32) {

    // TODO: update the code in this function
	double v = _v(P, opp_ratings, rating_list_size, opp_rds);
	P->vol = _newVol(P, opp_ratings, rating_list_size, opp_rds, scores, v);
	_preRatingRD(P);

	P->__rd = 1 / sqrt((1 / pow(P->__rd, 2)) + (1 / v));

	double tempSum = 0;
	for (int i = 0; i < rating_list_size; i++) {
		tempSum += _g(opp_rds[i]) * (scores[i] - _E(P, opp_ratings[i], opp_rds[i]));
	}
	P->__rating += pow(P->__rd, 2) * tempSum;
}


fn delta(player: &mut G2Player, opp_ratings: &f64, opp_rds: &f64,
                 scores: &f32) {

    // TODO: update the code in this function
}


/// Calculates the new volatility (σ') of the player with the rating deviation
/// `phi` (Φ), the volatility (σ) `sigma`, in a system with a system constant
/// (τ) `tau`, after a series of outcomes with a delta (∆) `delta` and a *v*
/// `v`. This function constitutes Step 5 of the Glicko-2 algorithm.
///
/// * `phi` the Glicko-2 scale value representing the player's rating deviation.
/// * `sigma` the value representing the player's volatility.
/// * `tau` the system constant.
/// * `v` the quantity *v* calculated in Step 3 from a series of outcomes
///     representing the estimated variance of the player.
/// * `delta` the quantity ∆ calculated in Step 4 from a series of outcomes
///     representing the estimated improvement in rating of the player.
/// Return:
/// * An f64 representing the player's new volatility σ'.
fn newVol(phi: f64, sigma: f64, tau: f64, v: f64, delta: f64) -> f64 {
    // Step 5, item 1
	let a: f64 = (sigma.powi(2)).ln();

    fn f(x: f64) -> f64 {
        let numerator1 = E.powf(x) * (delta.powi(2) - phi.powi(2) - v - E.powf(x));
        let denominator1 = 2.0 * (phi.powi(2) + v + E.powf(x)).powi(2);

        (numerator1 / denominator1) - ((x - a) / tau.powi(2))
    }
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
