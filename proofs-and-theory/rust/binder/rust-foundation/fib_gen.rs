use async_stream::stream;
use futures::Stream;
use futures::StreamExt;

async fn main() {
    let mut fib_stream = fibonacci().take(10);
    
    while let Some(value) = fib_stream.next().await {
        println!("{}", value);
    }
}

fn fibonacci() -> impl Stream<Item = u64> {
    stream! {
        let mut a = 0;
        let mut b = 1;

        loop {
            let next = a + b;
            a = b;
            b = next;
            yield a;
        }
    }
}