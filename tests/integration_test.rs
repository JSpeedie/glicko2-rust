use assert_float_eq::*;
use glicko2::*;

#[test]
fn glicko2_example() {
    let rating: f64 = 1500.0;
    let rd: f64 = 200.0;
    let sigma: f64 = 0.06;

    let tau: f64 = 0.5;

    let opp_ratings: [f64; 3] = [1400.0, 1550.0, 1700.0];
    let opp_rds: [f64; 3] = [30.0, 100.0, 300.0];
    let scores: [f64; 3] = [1.0, 0.0, 0.0];

    let mut mus: [f64; 3] = [0.0; 3];
    let mut phis: [f64; 3] = [0.0; 3];
    let mut gs: [f64; 3] = [0.0; 3];
    let mut es: [f64; 3] = [0.0; 3];

    // Check Step 2
    // Convert all player-of-interest rating and RD and all the opponent ratings and RDs to the
    // Glicko-2 scale.
    let mu = rating_to_g2_scale(rating);
    let phi = rd_to_g2_scale(rd);

    for i in 0..opp_ratings.len() {
        mus[i] = rating_to_g2_scale(opp_ratings[i]);
        phis[i] = rd_to_g2_scale(opp_rds[i]);
    }
    assert_f64_near!(mu, 0.0);
    assert_f64_near!(mus[0], -0.5756462492617337);
    assert_f64_near!(mus[1], 0.28782312463086684);
    assert_f64_near!(mus[2], 1.1512924985234674);
    assert_f64_near!(phi, 1.1512924985234674);
    assert_f64_near!(phis[0], 0.1726938747785201);
    assert_f64_near!(phis[1], 0.5756462492617337);
    assert_f64_near!(phis[2], 1.726938747785201);

    // Check `g()` calculations
    // Calculate g(Φ_j) for all phis
    for i in 0..gs.len() {
        gs[i] = _g(phis[i]);
    }
    assert_f64_near!(_g(phis[0]), 0.9954980064506083);
    assert_f64_near!(_g(phis[1]), 0.9531489778689763);
    assert_f64_near!(_g(phis[2]), 0.7242354780877525);

    // Check `E()` calculations
    // Calculate E(μ, μ_j, Φ_j) for all opponents
    for i in 0..es.len() {
        es[i] = _E(mu, mus[i], phis[i]);
    }
    assert_f64_near!(es[0], 0.6394677305521533);
    assert_f64_near!(es[1], 0.4318423561076679);
    assert_f64_near!(es[2], 0.3028407290952193);

    // Check *v* quantity calculations
    // Calculate v given (μ, μ_j, Φ_j)
    let v = _v(mu, &mus, &phis);
    assert_f64_near!(v, 1.7789770897239976);

    // Check delta (∆) quantity calculations
    // Calculate delta (∆) given (μ, μ_j, Φ_j, s_j, v)
    let delta = delta(mu, &mus, &phis, &scores, v);
    assert_f64_near!(delta, -0.4839332609836549);

    // Check the new volatility (σ') calculations
    // Calculate the new volatility (σ') given (Φ, σ, τ, v, ∆)
    let new_sigma = new_vol(phi, sigma, tau, v, delta);
    assert_f64_near!(new_sigma, 0.0599959842864885);

    // Check provisional phi (Φ*) calculations
    // Calculate the provisional phi (Φ*) given (Φ, σ')
    let prov_phi = update_rd(phi, new_sigma);
    assert_f64_near!(prov_phi, 1.1528546895801364);

    // Check new phi and mu (Φ', μ') calculations
    // Calculate the new phi and new mu (Φ', μ') given (μ, Φ*, μ_j, Φ_j, s_j, v)
    let (new_phi, new_mu) = update_player(mu, prov_phi, &mus, &phis, &scores, v);
    assert_f64_near!(new_phi, 0.8721991881307343);
    assert_f64_near!(new_mu, -0.20694096667525494);

    // Check the final result
    let new_rating = rating_to_original_scale(new_mu);
    let new_rd = rd_to_original_scale(new_phi);
    assert_f64_near!(new_rating, 1464.0506705393013);
    assert_f64_near!(new_rd, 151.51652412385727);
}
