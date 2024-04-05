import re

with open("src/bin/sokobond.rs") as f:
    code = f.read()


for test_case in re.finditer(r"test!\s*\{.*?\}", code, re.DOTALL):
    test_case = test_case.group(0)
    test_case = test_case.replace("\n", " ")

    test_case = test_case.split("{", 1)[1]
    test_case = test_case.split("}", 1)[0]
    test_case = test_case.split(",", 1)[1]

    _a, folder, _b, file, _c, test_case = test_case.split('"', 5)
    solutions = [
        solution for solution in re.findall(r"\"(.*?)\"", '"' + test_case) if solution
    ]

    print(folder, file, solutions)

    with open(f"data/sokobond/{folder}/{file}", "r") as f:
        content = f.read()

    content = content.strip()

    for solution in solutions:
        content += f"\n={solution}"
    content += "\n"

    print(content)

    with open(f"data/sokobond/{folder}/{file}", "w") as f:
        f.write(content)
