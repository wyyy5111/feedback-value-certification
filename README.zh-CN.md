# Feedback Value Certification

用于评估内生反馈下收缩型学习器有限时域期望价值的核心代码。

包含三个模块：

- 边界条件验证：合并相同投影增量，检验有理数词序列分离条件，计算有限状态块存活概率上界。
- 标量精确求值：前向构造可达包络，后向递推有理数价值分区，保留阈值点的实际标签。
- RPS 区间求值：在有理数猜拳与赢则保持、输则切换模型中，结合精确前缀、因果网格和定向整数舍入，输出确定性期望值包络。

## 安装

需要 Python 3.10 及以上版本。Python 核心使用标准库。C++ 模块需要 CMake 3.16 及以上版本，以及支持 `__int128` 的 C++17 编译器，例如 GCC 或 Clang。

```bash
git clone https://github.com/wyyy5111/feedback-value-certification.git
cd feedback-value-certification
python -m pip install .
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
```

## 使用

验证 RPS 投影分组的有限条件：

```bash
feedback-value certificate --input examples/rps_certificate.json
```

结果包含块长度 `1`、块存活概率上界 `179/180`，以及投影碰撞分组和分离数据。

计算标量模型的精确期望值：

```bash
feedback-value scalar --input examples/reactive_binary.json --policy tracking --horizon 6
```

结果为精确有理数 `7467729/9765625`。将策略改为 `contrarian` 可评估反向策略诱导的反馈规律。也可使用 `python -m feedback_value_certification` 调用相同接口。

计算 RPS 的确定性区间：

```bash
./build/rps-value --T 6 --G 8
```

结果以 JSON 输出，每个端点由整数分子和正分母表示。`BR_minus_mode` 给出两种期望值之差的区间。精确前缀长度为 `min(T, 8)`，后续网格递推在新的有限状态随机抽样之前选择共享子分片。整数舍入分母为 `2^40`。

## 输入与返回值

输入使用整数或 `"2/3"` 形式的有理数字符串。阈值点标签与开区间标签分别指定。响应和收益表保留当前收益与下一物理状态的联合依赖；标量模型的收益绝对值不超过一。

标量接口默认限制每个状态、每层最多 12,000 个切点，累计最多构造 200,000 个原子。达到限制时，Python 接口抛出 `EvaluationLimitError`，命令行返回 `RESOURCE_LIMIT` 和退出码 `2`；输入错误的退出码为 `1`。充分条件验证未通过会在结果中明确标记为无法确认该充分条件。

RPS 模块对应收缩因子 `49/50`、主选动作概率 `58/60`、首次索引优先的并列处理规则、学习器初始点 `(1/3,1/3,1/3)`，以及初始对手动作均匀分布的有理数模型。内存预检使用一 GiB 的保守分配上界。

## 基础检查

```bash
RPS_VALUE_EXECUTABLE="$PWD/build/rps-value" python -m unittest discover -s tests -v
```

检查覆盖短时域精确值、阈值标签、概率合法性、构造量限制、因果后续递推以及 RPS 端点回归。GitHub Actions 使用 Python 3.10 和 3.12 执行这些检查。

## 目录

```text
src/feedback_value_certification/  Python 核心算法与命令行接口
cpp/rps_value.cpp                 C++ RPS 区间算法
examples/                        两个最小输入示例
tests/                           单元与接口检查
```

项目使用 [MIT 许可证](LICENSE)。核心来源与模型归属见 [PROVENANCE.md](PROVENANCE.md)。固定版本源码可通过仓库的版本标签获取。
