test which suffix="":
    cargo build --release --bin {{which}}
    testit update --dry-run testit/{{which}}{{suffix}}.json

update which suffix="":
    cargo build --release --bin {{which}}
    testit update testit/{{which}}{{suffix}}.json