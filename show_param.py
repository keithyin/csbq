ALL_BASES = ["A", "C", "G", "T"]
import os


def decode_base_ctx(enc: int, ctx_len: int) -> str:
    stack = []

    while ctx_len > 0:
        v = enc & 0b11  # 取最低两位
        stack.append(ALL_BASES[v])
        enc >>= 2  # 右移两位
        ctx_len -= 1

    # 反转并返回字符串
    return "".join(reversed(stack))


def refcalled2prob(fname):
    refcalled2prob = {}

    with open(fname, "r") as file:
        for line in file:
            line = line.replace(" ", "").strip()  # 去掉空格和两端空白字符

            # 忽略注释行和不包含 '=' 的行
            if line.startswith("##") or "=" not in line:
                continue

            # 分割键和值
            left, right = line.split("=", 1)
            refcalled2prob[left] = float(right)

    return refcalled2prob


def main():
    ref2prob = refcalled2prob("model.param")
    results = []
    for k, prob in ref2prob.items():
        ctx = k.split("->", 1)[0]
        ctx = int(ctx)
        if k[-1] != ">":
            results.append(
                "{} ->{} {} = {}".format(
                    decode_base_ctx(ctx, ctx_len=3),
                    "+" if "+" in k else "",
                    decode_base_ctx(int(k[-1]), 1),
                    prob,
                )
            )
        else:
            results.append(
                "{} -> {} = {}".format(
                    decode_base_ctx(ctx, ctx_len=3),
                    " ",
                    prob,
                )
            )

    results = sorted(results)
    print("\n".join(results))


if __name__ == "__main__":
    main()
