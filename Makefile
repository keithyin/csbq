test:
	cargo test -- --nocapture

build:
	cargo build --release

install:
	cp target/release/csbq /usr/bin/

clean:
	rm -rf target