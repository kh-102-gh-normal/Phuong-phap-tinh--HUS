import math
import sympy as sp
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application,
                                        convert_xor)

X = sp.Symbol('x')
MAX_LAP = 10000
TF = standard_transformations + (implicit_multiplication_application, convert_xor)
LOCAL = {'x': X, 'e': sp.E, 'pi': sp.pi, 'ln': sp.log, 'tg': sp.tan, 'cotg': sp.cot,
         'arctg': sp.atan, 'arcsin': sp.asin, 'arccos': sp.acos}


class HamSo:
    def __init__(self, chuoi):
        self.expr = parse_expr(chuoi, transformations=TF, local_dict=LOCAL)
        self.d1_expr = sp.simplify(sp.diff(self.expr, X))
        self.d2_expr = sp.simplify(sp.diff(self.d1_expr, X))
        self.f = sp.lambdify(X, self.expr, 'math')
        self.d1 = sp.lambdify(X, self.d1_expr, 'math')
        self.d2 = sp.lambdify(X, self.d2_expr, 'math')


def an_toan(g, x):
    try:
        v = float(g(x))
        return v if math.isfinite(v) else None
    except Exception:
        return None


def dau(v):
    return (v > 0) - (v < 0)


def lay_mau(g, a, b, n=2000):
    return [an_toan(g, a + (b - a) * i / n) for i in range(n + 1)]


def chia_doi_tim(g, a, b, lan=60):
    ga = g(a)
    for _ in range(lan):
        c = (a + b) / 2
        gc = g(c)
        if ga * gc <= 0:
            b = c
        else:
            a, ga = c, gc
    return (a + b) / 2


def so_chu_so(eps):
    return max(6, int(math.ceil(-math.log10(eps))) + 3)


def fmt(v, k):
    return f"{v:.{k}f}"


def pstr(e):
    return str(e).replace("**", "^").replace("log", "ln")


def tieu_de(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)


def in_bang(tieu, hang):
    cot = len(tieu)
    rong = [max(len(str(tieu[j])), *(len(str(r[j])) for r in hang)) for j in range(cot)]
    duong = "+" + "+".join("-" * (w + 2) for w in rong) + "+"

    def dong(r):
        return "|" + "|".join(" " + str(r[j]).center(rong[j]) + " " for j in range(cot)) + "|"

    print(duong)
    print(dong(tieu))
    print(duong)
    for r in hang:
        print(dong(r))
    print(duong)


def lydo_khoang(hs, a, b):
    fa, fb = an_toan(hs.f, a), an_toan(hs.f, b)
    if fa is None or fb is None:
        return "hàm không xác định tại đầu mút"
    if fa * fb >= 0:
        return "f(a)·f(b) ≥ 0"
    v1 = lay_mau(hs.d1, a, b, 500)
    if any(v is None for v in v1):
        return "f'(x) không xác định trong khoảng"
    if not (all(v > 0 for v in v1) or all(v < 0 for v in v1)):
        return "f'(x) đổi dấu hoặc bằng 0 trong khoảng"
    v2 = lay_mau(hs.d2, a, b, 500)
    if any(v is None for v in v2):
        return "f''(x) không xác định trong khoảng"
    if not ((all(v >= -1e-12 for v in v2) or all(v <= 1e-12 for v in v2)) and any(abs(v) > 1e-12 for v in v2)):
        return "f''(x) đổi dấu trong khoảng"
    return None


def tim_khoang_ly_nghiem(hs, R=10):
    buoc = 0.01
    n = int(2 * R / buoc)
    xs = [round(-R + i * buoc, 10) for i in range(n + 1)]
    ys = [an_toan(hs.f, x) for x in xs]
    nghiem = []
    for i in range(n):
        if ys[i] is None or ys[i + 1] is None:
            continue
        if ys[i] == 0:
            nghiem.append(xs[i])
        elif ys[i] * ys[i + 1] < 0:
            r0 = chia_doi_tim(hs.f, xs[i], xs[i + 1])
            if abs(hs.f(r0)) < 1e-6:
                nghiem.append(r0)
    khoang, nhat_ky = [], []
    for r in nghiem:
        buoc_thu = []
        for s in (1, 0.1, 0.01, 0.001):
            a = round(math.floor(r / s) * s, 6)
            b = round(a + s, 6)
            ly_do = lydo_khoang(hs, a, b)
            buoc_thu.append((a, b, an_toan(hs.f, a), an_toan(hs.f, b), ly_do))
            if ly_do is None:
                if (a, b) not in khoang:
                    khoang.append((a, b))
                break
        nhat_ky.append((r, buoc_thu))
    return khoang, nhat_ky


def in_nhat_ky(nhat_ky):
    for i, (r, buoc_thu) in enumerate(nhat_ky, 1):
        print(f"\nNghiệm thứ {i} (nghiệm thô x ≈ {r:.4f}): tính f tại các đầu mút rồi thu hẹp khoảng")
        hang = []
        for a, b, fa, fb, ly_do in buoc_thu:
            fa_s = f"{fa:.4f}" if fa is not None else "?"
            fb_s = f"{fb:.4f}" if fb is not None else "?"
            kl = "THOẢ: f(a)f(b)<0, f', f'' không đổi dấu" if ly_do is None else "loại: " + ly_do
            hang.append([f"[{a:g}, {b:g}]", fa_s, fb_s, kl])
        in_bang(["khoảng", "f(a)", "f(b)", "kết luận"], hang)


def _gh(hs, diem, huong):
    try:
        L = sp.limit(hs.expr, X, diem, huong)
        t = str(sp.N(L, 6) if L.is_number and not L.has(sp.oo, sp.zoo) else L)
        return t.replace("zoo", "∞").replace("-oo", "-∞").replace("oo", "+∞")
    except Exception:
        return "?"


def _ngoai_mien(hs, x):
    return an_toan(hs.f, x) is None or an_toan(hs.d1, x) is None


def _bien_mien(hs, khong_ok, ok):
    for _ in range(50):
        m = (khong_ok + ok) / 2
        if _ngoai_mien(hs, m):
            khong_ok = m
        else:
            ok = m
    return ok


def bang_bien_thien(hs, R=10):
    buoc = 0.01
    n = int(2 * R / buoc)
    xs = [round(-R + i * buoc, 10) for i in range(n + 1)]
    ok = [not _ngoai_mien(hs, x) for x in xs]
    doan, i = [], 0
    while i <= n:
        if ok[i]:
            j = i
            while j < n and ok[j + 1]:
                j += 1
            doan.append((i, j))
            i = j + 1
        else:
            i += 1
    if not doan:
        print("Không tìm được miền xác định của hàm số trong [-10, 10].\n")
        return
    print("BẢNG BIẾN THIÊN")
    for (i0, i1) in doan:
        if i0 == 0:
            dau_doan = ("inf", -math.inf)
        else:
            b = _bien_mien(hs, xs[i0 - 1], xs[i0])
            dau_doan = ("bien", round(b, 6))
        if i1 == n:
            cuoi_doan = ("inf", math.inf)
        else:
            b = _bien_mien(hs, xs[i1 + 1], xs[i1])
            cuoi_doan = ("bien", round(b, 6))
        cuc_tri = []
        for k in range(i0, i1):
            d0, d1 = hs.d1(xs[k]), hs.d1(xs[k + 1])
            if abs(d0) < 1e-9:
                if k > i0 and hs.d1(xs[k - 1]) * d1 < 0:
                    cuc_tri.append(xs[k])
            elif d0 * d1 < 0:
                pt = chia_doi_tim(hs.d1, xs[k], xs[k + 1])
                if abs(hs.d1(pt)) < 1e-3:
                    cuc_tri.append(pt)
        diem = [dau_doan] + [("cuc", c) for c in cuc_tri] + [cuoi_doan]

        hx, hd, hy = [], [], []
        for idx, (loai, p) in enumerate(diem):
            if loai == "inf":
                hx.append("-∞" if p < 0 else "+∞"); hd.append("")
                hy.append(_gh(hs, -sp.oo if p < 0 else sp.oo, "+" if p < 0 else "-"))
            elif loai == "bien":
                hx.append(f"{p:g}"); hd.append("")
                hy.append(_gh(hs, sp.nsimplify(p), "+" if idx == 0 else "-"))
            else:
                hx.append(f"{p:.4f}"); hd.append("0"); hy.append(f"{hs.f(p):.4f}")
            if idx < len(diem) - 1:
                a, b = diem[idx][1], diem[idx + 1][1]
                if a == -math.inf and b == math.inf: m = 0.0
                elif a == -math.inf: m = b - 1
                elif b == math.inf: m = a + 1
                else: m = (a + b) / 2
                d = an_toan(hs.d1, m)
                sg = dau(d) if d is not None else 0
                hx.append(""); hd.append("+" if sg > 0 else "-" if sg < 0 else "?")
                hy.append("↗" if sg > 0 else "↘" if sg < 0 else "?")
        w = max(8, max(len(c) for c in hx + hd + hy) + 2)
        def dong(nhan, hang):
            return f"{nhan:>3} |" + "".join(c.center(w) for c in hang)
        if len(doan) > 1:
            print(f"  Trên khoảng ({dau_doan[1]:g}, {cuoi_doan[1]:g})".replace("inf", "∞"))
        print(dong("x", hx))
        print("-" * (5 + w * len(hx)))
        print(dong("y'", hd))
        print("-" * (5 + w * len(hx)))
        print(dong("y", hy))
        print()


def thong_so(hs, a, b):
    ga, gb = abs(hs.d1(a)), abs(hs.d1(b))
    return min(ga, gb), max(ga, gb)


def chieu_loi(hs, a, b):
    return "lồi (f'' > 0)" if hs.d2((a + b) / 2) > 0 else "lõm (f'' < 0)"


def ket_luan(ten, x, eps, a, b):
    d = max(1, int(math.ceil(-math.log10(eps))))
    print(f"\nKết luận ({ten}): nghiệm gần đúng trong [{a:g}, {b:g}] là")
    print(f"    x ≈ {x:.{d}f}   (giá trị tính được: {x:.{so_chu_so(eps)}f}, sai số < {eps:g})")


def chon_khoang(hs):
    tieu_de("1. KHẢO SÁT HÀM SỐ")
    print(f"f(x)   = {pstr(hs.expr)}")
    print(f"f'(x)  = {pstr(hs.d1_expr)}")
    print(f"f''(x) = {pstr(hs.d2_expr)}\n")
    bang_bien_thien(hs)
    tieu_de("2. KHOẢNG LY NGHIỆM")
    khoang, nhat_ky = tim_khoang_ly_nghiem(hs)
    if not khoang:
        print("Không tìm thấy khoảng ly nghiệm trong [-10, 10].")
        return None, None
    in_nhat_ky(nhat_ky)
    print("\nCác khoảng ly nghiệm tìm được:")
    for i, (a, b) in enumerate(khoang, 1):
        print(f"  {i}. [{a:g}, {b:g}]")
    i = input(f"\nChọn khoảng ly nghiệm (1-{len(khoang)}) [mặc định 1]: ").strip()
    return khoang[int(i) - 1 if i else 0]


def chia_doi(hs, a, b, eps):
    k = so_chu_so(eps)
    tieu_de("PHƯƠNG PHÁP CHIA ĐÔI")
    print("Bước 1: Khoảng ly nghiệm")
    print(f"    [a, b] = [{a:g}, {b:g}],  f(a) = {hs.f(a):.4f},  f(b) = {hs.f(b):.4f}  =>  f(a)·f(b) < 0")
    print("Bước 2: Công thức")
    print("    c = (a + b)/2;  nếu f(a)·f(c) < 0 thì b := c, ngược lại a := c")
    print("    Sai số sau n bước: |x_n - x*| ≤ (b - a)/2^n")
    n_du = max(1, math.ceil(math.log2((b - a) / eps)))
    print(f"    Cần (b - a)/2^n < {eps:g}  =>  n ≥ {n_du}")
    print("Bước 3: Thực hiện chia đôi")
    hang, n = [], 0
    while n < MAX_LAP:
        n += 1
        c = (a + b) / 2
        fa, fc = hs.f(a), hs.f(c)
        sai_so = (b - a) / 2
        hang.append([n, fmt(a, k), fmt(b, k), fmt(c, k), f"{fc:.3e}", f"{sai_so:.3e}"])
        if fc == 0 or sai_so < eps:
            break
        if fa * fc < 0:
            b = c
        else:
            a = c
    in_bang(["n", "a", "b", "c=(a+b)/2", "f(c)", "sai số (b-a)/2"], hang)
    return c, hang


def main():
    s = input("Nhập hàm số f(x) (vd: x^3 - 3x + 1): ").strip()
    hs = HamSo(s)
    eps = float(input("Nhập sai số yêu cầu (vd: 1e-4): "))
    a, b = chon_khoang(hs)
    if a is None:
        return
    x, _ = chia_doi(hs, a, b, eps)
    if x is not None:
        ket_luan("Chia đôi", x, eps, a, b)


if __name__ == "__main__":
    main()
