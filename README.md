# GNSS RTK Relative Positioning

一个基于 C++ 的 GNSS 相对定位实验项目，包含 SPP/RTK 相关计算流程，支持：

- 文件离线模式（PPK）
- 网络 socket 实时读流模式（RTK）

当前工程为 Visual Studio 解决方案，入口位于 `project/project.cpp`。

## 项目结构

- `project.sln`：Visual Studio 解决方案
- `project/`：核心源代码目录
  - `project.cpp`：程序入口，切换离线/网络模式
  - `RTK*.cpp/.h`：RTK 读取、建模与求解
  - `Spp.cpp`：单点定位相关流程
  - `decoding.cpp/.h`：原始数据解码
  - `matrix.cpp/.h`、`lambda.cpp/.h`：矩阵与整数模糊度相关计算
  - `sockets.cpp/.h`：Windows Socket 通信

## 开发环境

- Windows 10/11
- Visual Studio 2022（MSVC v143 工具集）
- C++ 控制台应用（已在 `project.vcxproj` 配置）
- 依赖系统库：`WS2_32.lib`

## 编译与运行

1. 使用 Visual Studio 打开 `project.sln`
2. 选择 `Debug|x64`（或其他你需要的配置）
3. 生成并运行项目

## 运行模式说明

`project/project.cpp` 中预留了两种模式：

### 1) 网络实时模式（默认）

当前代码默认调用：

- `RTKSockets(Sect1, Sect2);`

并使用两个测站的 IP/端口配置（在 `main()` 中可修改）：

- Base：`47.114.134.129:7190`
- Rover：`8.140.46.126:3002`

### 2) 文件离线模式（PPK）

在 `main()` 中可使用文件模式（示例代码已注释）：

- 打开基准站与流动站二进制观测文件
- 调用 `PPKFile(...)`

如需启用离线模式：

1. 取消 `FILE*` 和 `PPKFile(...)` 对应注释
2. 注释掉 `RTKSockets(...)`
3. 确认数据文件路径正确

## 数据文件与输出

仓库中包含部分示例二进制数据（`.bin`）。  
运行过程中可能会生成/使用如下结果文件（按代码逻辑而定）：

- `RTKrecord.txt`
- `RTKFloat.txt`
- `RTKRes.csv`
- `BaseSPPResult.csv`
- `RoverSPPResult.csv`

