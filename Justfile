test which suffix="":
    cargo build --release --bin {{which}}
    testit update testit/{{which}}{{suffix}}.json
