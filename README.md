# glicko2-rust

An implementation of the Glicko-2 algorithm written in Rust.

## Using the Library

If you have a Rust project in which you would like to make use of
this library, here's how you can set it up locally.

### 1. Clone this library

```bash
# If you want to push to the repo
git clone git@github.com:JSpeedie/glicko2-rust.git glicko2-rust
# If you just want a copy of the code
git clone https://www.github.com/JSpeedie/glicko2-rust glicko2-rust
```

### 2. Add the library as a dependency in your project

Edit the dependency section of the `Cargo.toml` file in your project
to contain the following line:

```toml
[dependencies]
glicko2 = { path = "../glicko2-rust/" } # Depend on the library as it exists in the local filesystem
```

### 3. Add a `use` statement in your project

Assuming you want to make use of this library in the `main.rs` file of
your project, you could add this near the top of your `main.rs` file:

```rust
use glicko2::*;
```

Make sure the module you use in this line of code (i.e. `use glicko2 ...`)
matches the dependency name in your Cargo.toml (`glicko2 = ...`). The name
does not need to match the name of the local directory.

Now that you are able to make use of the library, you're ready to make some
Glicko-2 calculations! An example is provided in the documentation that walks
you through an example series of function calls representing the calculation
of a player's new rating and rating deviation given a series of games. Take
a look at the section below on library documentation to learn more.

&nbsp;
&nbsp;

## Testing the Library

This library has both unit tests and integration tests. You can run them as you
would for any Rust project using `cargo test`. Here's what the output may look
like:

```bash
cargo test
```
```
    Finished `test` profile [unoptimized + debuginfo] target(s) in 0.00s
     Running unittests src/lib.rs (target/debug/deps/glicko2-61df8c83d8ed253e)

running 11 tests
test unit_tests::e_function ... ok
test unit_tests::calculate_new_rating_and_rd ... ok
test unit_tests::new_volatility ... ok
test unit_tests::g_function ... ok
test unit_tests::new_pre_rating_period_value ... ok
test unit_tests::scale_rating_to_original_scale ... ok
test unit_tests::delta_quantity ... ok
test unit_tests::scale_rd_to_g2_scale ... ok
test unit_tests::scale_rd_to_original_scale ... ok
test unit_tests::v_quantity ... ok
test unit_tests::scale_rating_to_g2_scale ... ok

test result: ok. 11 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

     Running tests/integration_test.rs (target/debug/deps/integration_test-f2eba1d7ef98fb33)

running 1 test
test integration_tests::glicko2_example ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

   Doc-tests glicko2

running 0 tests

test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s

```

&nbsp;
&nbsp;

## Library Documentation

This library is complete with documentation. You can view the documentation as
you would for any Rust project:

```bash
cargo doc --open
```

The root HTML page for the documentation should open up in your default web
browser.

&nbsp;
&nbsp;

## Linting

This project is linted by `clippy`. The process for linting is quite simple.
First, we want to make sure `clippy` is installed:

```bash
rustup component list | grep clippy
```
```
clippy-x86_64-unknown-linux-gnu (installed)
```

If the output of this command is empty, we need to install clippy. We can do so
with the following commands:

```bash
rustup update
rustup component add clippy
```

Once we have clippy installed, we can lint our program with the following
command:

```bash
cargo clippy
```
```
    Checking glicko2 v0.1.0 (/home/me/glicko2-rust)
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 0.07s
```

If the code passed all lints, the output will look similar to the example
output shown above.

To check the project more thoroughly, you can use:

```bash
cargo clippy -- -W clippy:pedantic
```

This typically lints a little too aggressively for my taste, but it's worth
checking every once in a while.
