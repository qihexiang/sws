pub mod definitions;
mod tokenizer;

fn create_counter() -> impl FnMut() -> usize {
    let mut count = 0usize;
    let counter = move || {
        count += 1;
        count
    };
    counter
}

fn adder(a: usize) -> impl Fn(usize) -> usize {
    move |b| a + b
}

#[test]
fn test() {
    let mut counter = create_counter();
    for _ in 0..5 {
        println!("{}", counter());
    }
    println!("{}", adder(1)(2))
}
