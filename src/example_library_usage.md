# Example Library Usage

Let's follow the example given in the Glicko-2 pdf. We have a player
(rating=1500, rd=200, volatility=0.06) who plays 3 games against 3 opponents:

* opp_1 (rating=1400, rd=30)
* opp_2 (rating=1550, rd=100)
* opp_3 (rating=1700, rd=300)

The first game our player-of-interest wins. The games against opponent 2 and 3,
our player-of-interest loses. For this example, let's set the system constant τ
at 0.5. We can represent these elements of the scenario in code like so:

```text
// The representation of our player-of-interest
let rating: f64 = 1500.0;
let rd: f64 = 200.0;
let sigma: f64 = 0.06;

// The system constant
let tau: f64 = 0.5;

// Arrays represent the ratings and rating deviations of the opponents our player-of-interest
// played, as well as the score of the games.
let opp_ratings: [f64; 3] = [1400.0, 1550.0, 1700.0];
let opp_rds: [f64; 3] = [30.0, 100.0, 300.0];
let scores: [f64; 3] = [1.0, 0.0, 0.0];
```

Next we need to convert from the original scale to the Glicko-2 scale the
rating and rating deviation numbers for both the player-of-interest as well as
all the opponents .

```text
let mut mus: [f64; 3] = [0.0; 3];
let mut phis: [f64; 3] = [0.0; 3];

// Convert the rating and rating deviation of our player-of-interest
let mu = rating_to_g2_scale(rating);
let phi = rd_to_g2_scale(rd);

// Convert the ratings and rating deviations of all the opponents
for i in 0..opp_ratings.len() {
    mus[i] = rating_to_g2_scale(opp_ratings[i]);
    phis[i] = rd_to_g2_scale(opp_rds[i]);
}
```

Next we need to get the output from two helper functions in the Glicko-2
algorithm: the *g()* function and the *E()* function. We need to get the output
from these functions once for each opponent our player-of-interest played. The
code for doing this could look like this:

```text
let mut gs: [f64; 3] = [0.0; 3];
let mut es: [f64; 3] = [0.0; 3];

// Calculate g(Φ_j) for all opponent phis
for i in 0..gs.len() {
    gs[i] = _g(phis[i]);
}

// Calculate E(μ, μ_j, Φ_j) for all opponents
for i in 0..es.len() {
    es[i] = _E(mu, mus[i], phis[i]);
}
```

Next we need to calculate the *v* quantity for this series of outcomes. We can
do that like so:

```text
// Calculate *v* given (μ, μ_j, Φ_j)
let v = _v(mu, &mus, &phis);
```

Once we have the *v* quantity, we need to calculate the delta (∆) quantity:

```text
// Calculate delta (∆) given (μ, μ_j, Φ_j, s_j, v)
let delta = delta(mu, &mus, &phis, &scores, v);
```

After that we calculate the new volatility for the player-of-interest:

```text
// Calculate the new volatility (σ') given (Φ, σ, τ, v, ∆)
let new_sigma = new_vol(phi, sigma, tau, v, delta);
```

Then we calculate the new pre-rating period value (Φ*) for the
player-of-interest:

```text
// Calculate the pre-rating period value (Φ*) given (Φ, σ')
let phi_star = pre_rating_period_value(phi, new_sigma);
```

The second last step is to calculate the new rating and rating deviation for our
player-of-interest:


```text
// Calculate the new phi and new mu (Φ', μ') given (μ, Φ*, μ_j, Φ_j, s_j, v)
let (new_phi, new_mu) = new_rating_and_rd(mu, phi_star, &mus, &phis, &scores, v);
```

And the final step is to convert that new rating and rating deviation to the
original scale:

```text
let new_rating = rating_to_original_scale(new_mu);
let new_rd = rd_to_original_scale(new_phi);
```

All done, our player-of-interest's new Glicko-2 data is (rating=1464.0508,
rd=151.5165, volatility=0.059996) as captured by the `new_rating`, `new_rd`,
and `new_sigma` variables.
