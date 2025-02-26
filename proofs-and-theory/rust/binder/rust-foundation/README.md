# Rust Foundation Guides

## 1. Closures
* Rust closures are anonymous functions you can save in a variable or pass as arguments to other functions.
* You can create the closure in one place, and call the closure elsewhere to evaluate in a different context
* (!) Unlike functions, closures can capture values from the scope in which they're defined.

## 2. Generator Functions
* Generator functions are a type of function that, instead of returning a single value, can yield multiple values over time. In Rust, genfuncs are implemented using coroutines, which are a generalization of subroutines that can be paused and resumed multiple times during their execution.

* A generator function can be thought of as a state machine that encapsulates the state of an iterator. Once a generator is paused, it can be resumed from the same state it was paused in, allowing it to continue producing values from where it left off.

* The example in this repo uses async-stream and fibonacci as a simple example