#![allow(unused)]

fn main() {
    use std::collections::HashMap;

    // Type inference lets us omit an explicit type signature (which
    // would be `HashMap<String, String>` in this example).
    let mut book_reviews = HashMap::new();

    // Review some books.
    book_reviews.insert(28, "My favorite book.".to_string());

    book_reviews.insert(888, "Hugo Ledoux".to_string());

    let mykey: usize = 28;
    if book_reviews.contains_key(&mykey) {
        println!("yeah");
    }

    book_reviews.remove(&29);

    // Look up the value for a key (will panic if the key is not found).
    // println!("Review for Jane: {}", book_reviews["Pride and Prejudice"]);

    // Iterate over everything.
    for (book, review) in &book_reviews {
        println!("{}: \"{}\"", book, review);
    }
}
