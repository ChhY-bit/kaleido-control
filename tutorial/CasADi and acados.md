# CasADi/acados 安装指南（MATLAB）

- **什么是 `CasADi` ？**
    符号建模、自动微分、代码生成及调用优化求解器的数值优化工具框架。

- **什么是 `acados` ？**
    利用最优控制问题结构，高效求解 MPC 等问题的框架，通常通过 CasADi 定义模型并生成所需函数。

## 安装和使用 CasADi

官方网站： https://web.casadi.org/get/

### 获取 `.zip` 安装包

从 GitHub 仓库 (https://github.com/casadi/casadi/releases) 获取最新版本，例如：

- [Linux 版本 `casadi-3.8.1-linux64-matlab2018b.zip`](https://github.com/casadi/casadi/releases/download/3.8.1/casadi-3.8.1-linux64-matlab2018b.zip) 
- [Win 版本 `casadi-3.8.1-windows64-matlab2018b.zip`](https://github.com/casadi/casadi/releases/download/3.8.1/casadi-3.8.1-windows64-matlab2018b.zip)

### 建议安装路径

下载完成后解压，并复制/移动到合适的路径，建议：

- 进入 MATLAB 的 Toolbox 路径，通常是 `your_matlab/toolbox/`
- 新建文件夹 `CasADi/`，将解压后的文件夹复制/移动至其中
- 最终路径大致为 `your_matlab/toolbox/CasADi/casadi-3.8.1-windows64-matlab2018b`

### 执行安装

启动MATLAB，将工作目录切换为上述路径（能看到有名为`+casadi`的子文件夹），运行：

```matlab
addpath(pwd)
savepath
```

即可完成安装。

可用以下片段测试是否正常工作：

```matlab
import casadi.*
x = MX.sym('x')
disp(jacobian(sin(x),x))
```

### 快速使用 CasADi

> 此处展示 CasADi 在 **二次规划(QP)** 和 **非线性规划(NLP)** 中的使用方法。