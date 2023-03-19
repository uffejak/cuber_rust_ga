# cuber_rust_ga
Simple genetic algorithm for Cuber project. 

Contains ideas borrowed from gradient descent and simulated annealing.


# Conditional compilation

Added conditional defines/features to Cargo TOML:
[features]
model_equvalent_circuit_2nd_order = []
model_equvalent_circuit_2nd_order_with_decay = []

"# features by default"
default = ["model_equvalent_circuit_2nd_order"]

“main.rs” will then include different sections of code based on features you choose.

So when you need to implement new Rust models search for sections such as:
"#[cfg(feature = "model_equvalent_circuit_2nd_order_with_decay")]"
To see where you need to insert your code/ copy+modify existing code

See also:
https://web.mit.edu/rust-lang_v1.25/arch/amd64_ubuntu1404/share/doc/rust/html/book/first-edition/conditional-compilation.html


# Funding
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under Grant Agreement no. 875605.
The results presented on this website reflect only CuBER’s view. The European Commission is not responsible for
any use that may be made of the information it contains.

# License
This project uses The boost license.
https://www.boost.org/users/license.html
