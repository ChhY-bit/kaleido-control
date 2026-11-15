# DeePC MATLAB-APIs

## APIs

### Methods

#### deepc(N,T_ini)
- `N` - 预测步数
- `T_ini` - 初始长度
- 作用：初始化 DeePC 控制器
- 示例：
    ```matlab
    % 用20步历史数据向前预测30步
    N = 20;
    T_ini = 30;
    controller = ctr.deepc(N,T_ini);
    ```

#### deepc.setParams(Q,R,lambda_u,lambda_y,lambda_a)
- `Q` - 跟踪代价权重
- `R` - 控制代价权重
- `lambda_u` - 输入松弛项正则化系数（默认 `1e6`）
- `lambda_y` - 输出松弛项正则化系数（默认 `1e4`）
- `lambda_a` - alpha项松弛正则化系数（默认 `1e2`）
- 作用：设置控制器参数
- 示例：
    ```matlab
    controller.setParams(Q,R,lambda_u,lambda_y,lambda_a);
    ```

#### deepc.init(u_d,y_d)
- `u_d` - 输入数据库
- `y_d` - 输出数据库
- 作用：初始化数据库，完成离线计算
- 示例：
    ```matlab
    % 数据获取部分略
    % u_d = ...
    % y_d = ...
    controller.init(u_d,y_d);  % 传入数据
    ```

#### (optional) deepc.samplify(ratio,ifplot)
- `ratio` - 奇异值截断比（默认 `10`）
- `ifplot` - 是否绘制 SVD 结果
- 作用：对数据信息进行降维，剔除噪声，保留最少的生成向量
- 注意：必须在 `.init()` 之后使用
- 示例：
    ```matlab
    controller.samplify(10,flase);  % 不开启绘图
    ```

#### deepc.setup(solver_name,settings)
- `solver_name` - 使用的求解器名称
- `settings` - 求解器的设置属性
- 作用：设定所采用的求解器，求解预热
- 注意：必须在 `.init()` 或 `.samplify()` 之后使用
- 示例：
    ```matlab
    controller.init(u_d,y_d);
    controller.samplify(10,flase);
    controller.setup('osqp',[]);    % 使用OSQP求解器，一般 [] 表示默认设定
    ```

#### deepc.update(u_ini,y_ini,u_r,y_r,Y_lb,Y_ub,U_lb,U_ub)
- `u_ini` - 长度为 `T_ini` 的初始输入数据
- `y_ini` - 长度为 `T_ini` 的初始输出数据
- `u_r` - 长度为 `N` 的目标控制量
- `y_r` - 长度为 `N` 的目标输出量
- `Y_lb` - 输出约束下界
- `Y_ub` - 输出约束上界
- `U_lb` - 控制约束下界
- `U_ub` - 控制约束上界
- 作用：在控制循环中实时反馈更新
- 示例：见 `.solve()`

#### deepc.solve()
- 作用：求解最优控制量并计算最优代价
- 示例：
    ```matlab
    for i = 1:1000      % 以1000次控制为例
        % 其他代码略
        controller.update(u_ini,y_ini,u_r,y_r,Y_lb,Y_ub,U_lb,U_ub)
        controller.solve()
    end
    ```

### Properties

## Example