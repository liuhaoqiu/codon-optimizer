# AGENTS.md

本文件用于指导在本仓库内工作的智能编码代理。
内容基于当前仓库结构（核心脚本为 `check_cds.py`）。

## 1）仓库概览

- 主要语言：Python 3
- 主入口脚本：`check_cds.py`
- 默认输入数据：`cds_from_genomic.fna`
- 主要输出产物：`cds_codon_report.html`
- 构建系统：未检测到
- Lint 配置：未检测到
- 测试配置：未检测到

## 2）Cursor / Copilot 规则检查

已检查以下文件/目录，当前均不存在：

- `.cursorrules`
- `.cursor/rules/`
- `.github/copilot-instructions.md`

若后续新增上述规则文件，优先级应高于本文件。

## 3）环境准备

建议所有本地运行与代理执行都使用虚拟环境。

### Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
```

### bash

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

### 运行时依赖

脚本依赖如下第三方包：

- `matplotlib`
- `numpy`
- `pandas`

当前脚本内含“缺包自动安装”逻辑。
为保证可复现性（CI 或代理批量执行），建议显式安装：

```bash
python -m pip install matplotlib numpy pandas
```

## 4）Build / Lint / Test 命令

当前仓库没有正式的构建流水线。

### 运行主程序（等价于集成运行）

```bash
python check_cds.py
```

预期结果：

- 读取并解析 `cds_from_genomic.fna`
- 在终端打印校验统计摘要
- 生成 `cds_codon_report.html`

### Build 命令现状

- 目前不存在独立 build 步骤。
- 可将 `python check_cds.py` 视为端到端执行命令。

### Lint 命令（推荐临时使用）

```bash
python -m pip install ruff
ruff check check_cds.py
```

### 格式检查命令（推荐临时使用）

```bash
python -m pip install black
black --check check_cds.py
```

### 测试命令现状

- 当前仓库尚无自动化测试文件。
- 若后续引入 `pytest`，可执行全部测试：

```bash
pytest
```

### 单测单条运行（重点）

当存在 pytest 用例时，运行单条测试可用：

```bash
pytest tests/test_file.py::test_name
pytest tests/test_file.py -k "test_name"
pytest -k "substring"
```

### 当前仓库可执行的轻量校验

```bash
python -m py_compile check_cds.py
python check_cds.py
```

## 5）代码风格与工程规范

若仓库后续新增格式化/Lint 配置文件，以配置为准。
在此之前，遵循 `check_cds.py` 的现有风格。

### Python 版本与兼容性

- 目标语义：Python 3.9+
- 运行环境：标准 CPython

### Import 规范

- 按顺序分组：标准库、第三方、本地模块。
- 保持导入顺序稳定，避免无意义改动。
- 非必要不要重复导入。
- 若因启动顺序必须延迟导入，需在代码中说明原因。

### 格式化规范

- 使用 4 空格缩进，禁止 Tab。
- 建议行宽约 88-100 字符。
- 逻辑块之间保留一个空行。
- 多行字面量可使用尾随逗号，提升 diff 可读性。

### 类型注解规范

- 对新增/修改的公共函数尽量补充类型注解。
- 可复用辅助函数建议显式返回类型。
- 优先使用内建泛型（如 `list[str]`、`dict[str, int]`）。

### 命名规范

- 函数/变量：`snake_case`
- 类名：`PascalCase`
- 常量：`UPPER_SNAKE_CASE`
- 生物信息学相关变量命名应语义清晰（如 `valid_stop_codons`）。

### 字符串与模板

- 字符串插值优先使用 f-string。
- 大段 HTML/CSS/JS 建议使用三引号块。
- 写入 HTML 前，必须对动态数据进行转义。
- 面向用户的 HTML 文案统一使用英文，不使用中文界面文本。

### 错误处理规范

- 尽量捕获具体异常类型。
- 避免裸 `except:`。
- 错误信息应包含上下文且可执行（可定位、可修复）。
- 对关键输入缺失应尽早失败（fail fast）。

### 校验与数据处理

- 序列应统一 `upper()` 并去除空白字符。
- 明确校验规则：长度是否为 3 倍数、起始/终止密码子、ATCG 字符集。
- 校验输出尽量结构化，保证可重复与可比对。

### I/O 与副作用

- 文本读写显式指定 UTF-8。
- 默认路径应可预期，避免隐式行为。
- 尽量分离解析、校验、统计、报告渲染四类职责。

### 注释与文档字符串

- Docstring 简洁、以行为说明为主。
- 仅对不直观逻辑补充注释。
- 修改代码时同步清理过期注释。

### 依赖管理

- 建议后续增加显式依赖文件（如 `requirements.txt` 或 `pyproject.toml`）。
- 若保留自动安装逻辑，失败提示必须清晰并含修复建议。

### 性能与可扩展性

- 当前 FASTA 解析为整文件读入内存。
- 若数据规模增大，重构时优先考虑流式/分块处理。

## 6）代理工作流程要求

- 默认进行最小、聚焦修改，除非明确要求重构。
- 非需求驱动情况下，不主动改变既有行为。
- 代码变更后至少执行语法检查与主流程运行验证。
- 汇报需包含：改了什么、为什么改、如何验证。
