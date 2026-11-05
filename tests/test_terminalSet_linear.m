function tests = test_terminalSet_linear
%TEST_TERMINALSET_LINEAR 验证夹逼关系、闭环不变性和公共接口。
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    root = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.originalPath = path;
    addpath(fullfile(root,'lib'));
end

function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end

function testRotatedEllipsoid(testCase)
    M = [2 .6;.6 1];
    R = chol(M);
    theta = .7;
    rotation = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    At = R\(.88*rotation)*R;
    verifyGeometry(testCase,At,M,2.3);
end

function testThreeDimensions(testCase)
    M = [3 .4 -.2;.4 2 .3;-.2 .3 1];
    R = chol(M);
    At = R\[.6 -.2 0;.2 .6 0;0 0 .72]*R;
    verifyGeometry(testCase,At,M,.4);
end

function testWeakContraction(testCase)
    angle = .4;
    At = .98*[cos(angle) -sin(angle);sin(angle) cos(angle)];
    verifyGeometry(testCase,At,eye(2),1);
end

function testDefaultBetaMatchesExplicitMinimum(testCase)
    At = [.65 .15;0 .7];
    M = [2 .3;.3 1];
    R = chol(M);
    beta = norm(R*At/R,2);
    [lb1,ub1,E1] = ctr.terminalSet_linear(At,M,1);
    [lb2,ub2,E2] = ctr.terminalSet_linear(At,M,1,beta);
    verifyEqual(testCase,E2,E1,'AbsTol',1e-12);
    verifyEqual(testCase,lb2,lb1,'AbsTol',1e-12);
    verifyEqual(testCase,ub2,ub1,'AbsTol',1e-12);
end

function testCustomBetaControlsInnerEllipsoid(testCase)
    % 实际最小收缩因子为 .5；显式要求包含 .95 倍椭球。
    verifyGeometry(testCase,.5*eye(2),eye(2),2,.95);
end

function testInvalidBetaReportsInterval(testCase)
    invalid = {.79,1,1.1,-.1,NaN,Inf,-Inf,[],[.8 .9],.9+1i,'0.9',true,struct()};
    lower = sprintf('[%.17g, 1)',norm(.8*eye(2),2));
    for i = 1:numel(invalid)
        caught = [];
        try
            ctr.terminalSet_linear(.8*eye(2),eye(2),1,invalid{i});
        catch exception
            caught = exception;
        end
        assertNotEmpty(testCase,caught);
        verifyEqual(testCase,caught.identifier,'terminalSet_linear:InvalidBeta');
        verifySubstring(testCase,caught.message,lower);
    end
end

function testBetaValidatedBeforeSpecialCases(testCase)
    verifyError(testCase,@() ctr.terminalSet_linear(.8,1,1,.7), ...
        'terminalSet_linear:InvalidBeta');
    verifyError(testCase,@() ctr.terminalSet_linear(.8*eye(2),eye(2),0,1), ...
        'terminalSet_linear:InvalidBeta');
    verifyGeometry(testCase,zeros(2),eye(2),1,0);
end

function testSinglePrecisionBetaBoundary(testCase)
    % single(.7) 的 double 值略小于 .7，不能因混合精度比较而误接受。
    verifyError(testCase,@() ctr.terminalSet_linear(.7*eye(2),eye(2),1,single(.7)), ...
        'terminalSet_linear:InvalidBeta');
    verifyGeometry(testCase,.5*eye(2),eye(2),1,single(.9));
end

function testFourDimensionalCoplanarFacets(testCase)
    [lb,ub,E] = ctr.terminalSet_linear(.6*eye(4),eye(4),1);
    verifyEqual(testCase,size(E,2),4);
    verifyEqual(testCase,lb,-ub);
    verifyLessThanOrEqual(testCase,.6*vecnorm(E,2,2),ub+1e-9);
    % 四维凸包的共面面片回归：用独立 LP 检查坐标轴和各卦限方向。
    directions = [eye(4);-eye(4);2*(dec2bin(0:15,4)-'0')-1];
    options = optimoptions('linprog','Display','none');
    for i = 1:size(directions,1)
        [point,~,flag] = linprog(-directions(i,:)',[E;-E],[ub;-lb], ...
            [],[],[],[],options);
        verifyGreaterThan(testCase,flag,0);
        verifyLessThanOrEqual(testCase,point'*point,1+1e-8);
        verifyLessThanOrEqual(testCase,abs(E*(.6*point)),ub+1e-8);
    end
end

function testSingularClosedLoop(testCase)
    verifyGeometry(testCase,[.8 .1;0 0],[2 .2;.2 1],1);
end

function testZeroDynamics(testCase)
    verifyGeometry(testCase,zeros(2),[1 .2;.2 2],3);
end

function testScaling(testCase)
    At = [.65 .15;0 .7];
    M = [2 .3;.3 1];
    [lb1,ub1,E1] = ctr.terminalSet_linear(At,M,1);
    [lb2,ub2,E2] = ctr.terminalSet_linear(At,M,9);
    verifyEqual(testCase,E2,E1,'AbsTol',1e-12);
    verifyEqual(testCase,lb2,3*lb1,'AbsTol',1e-12);
    verifyEqual(testCase,ub2,3*ub1,'AbsTol',1e-12);
end

function testOneDimension(testCase)
    [lb,ub,E] = ctr.terminalSet_linear(-.9,4,9);
    verifyEqual(testCase,E,2);
    verifyEqual(testCase,lb,-3);
    verifyEqual(testCase,ub,3);
end

function testZeroLevel(testCase)
    [lb,ub,E] = ctr.terminalSet_linear(.8*eye(3),eye(3),0);
    verifyEqual(testCase,E,eye(3));
    verifyEqual(testCase,lb,zeros(3,1));
    verifyEqual(testCase,ub,zeros(3,1));
end

function testRejectsInvalidLyapunovMatrix(testCase)
    verifyError(testCase,@() ctr.terminalSet_linear(.5*eye(2), ...
        [1 .2;0 1],1),'terminalSet_linear:NotSymmetric');
    verifyError(testCase,@() ctr.terminalSet_linear(.5*eye(2), ...
        diag([1 -1]),1),'terminalSet_linear:NotPositiveDefinite');
    % Schur 稳定本身不足以使给定的 M 成为一步收缩椭球。
    verifyError(testCase,@() ctr.terminalSet_linear([.5 2;0 .5], ...
        eye(2),1),'terminalSet_linear:NotContractive');
    verifyError(testCase,@() ctr.terminalSet_linear(eye(2), ...
        eye(2),1),'terminalSet_linear:NotContractive');
end

function testInvalidInputShapes(testCase)
    verifyError(testCase,@() ctr.terminalSet_linear(ones(2,3),eye(2),1), ...
        'MATLAB:terminalSet_linear:expectedSquare');
    verifyError(testCase,@() ctr.terminalSet_linear(.5*eye(2),eye(2),-1), ...
        'MATLAB:terminalSet_linear:expectedNonnegative');
end

function verifyGeometry(testCase,At,M,alpha,varargin)
    [lb,ub,E] = ctr.terminalSet_linear(At,M,alpha,varargin{:});
    n = size(At,1);
    verifySize(testCase,E,[numel(ub) n]);
    verifySize(testCase,ub,[size(E,1) 1]);
    verifyEqual(testCase,lb,-ub);
    verifyTrue(testCase,all(isfinite(E(:))) && all(isfinite(ub)));
    verifyGreaterThan(testCase,ub,zeros(size(ub)));

    % beta*E0 包含于 P：椭球的线性支撑函数给出逐行的精确最大值。
    R = chol(M);
    beta = norm(R*At/R,2);
    if ~isempty(varargin)
        beta = double(varargin{1});
    end
    support = beta*sqrt(alpha)*vecnorm(E/R,2,2);
    verifyLessThanOrEqual(testCase,support,ub+1e-9);

    % 独立枚举半空间交点，而非复用实现中的凸包/极集认证代码。
    A = [E;-E];
    b = [ub;-lb];
    combinations = nchoosek(1:size(A,1),n);
    vertices = zeros(0,n);
    for i = 1:size(combinations,1)
        active = A(combinations(i,:),:);
        if rcond(active) < 1e-11
            continue
        end
        point = active\b(combinations(i,:));
        if all(A*point <= b+1e-8*max(1,norm(point)))
            vertices(end+1,:) = point'; %#ok<AGROW>
        end
    end
    verifyNotEmpty(testCase,vertices);
    verifyEqual(testCase,rank(vertices),n);
    values = sum((vertices*M).*vertices,2);
    verifyLessThanOrEqual(testCase,values,alpha*(1+1e-8)*ones(size(values)));
    next = E*At*vertices';
    verifyLessThanOrEqual(testCase,max(next,[],2),ub+1e-8);
    verifyGreaterThanOrEqual(testCase,min(next,[],2),lb-1e-8);
end
