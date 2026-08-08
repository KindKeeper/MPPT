%% run_all.m — MPPT 仓库一键运行脚本
% 功能：依次运行 代码/ 目录下所有可独立运行的 MPPT 算法脚本，验证仓库内容可复现。
% 运行环境：MATLAB 或 Octave（无图形界面环境下可自动跳过绘图）
% 用法：在仓库根目录执行  run_all  即可。
%
% 说明：每个脚本以独立子进程方式运行，互不干扰，
%       脚本内部的 clear/clc/figure 不会影响其他脚本或本主脚本。

% 获取本脚本所在目录（确保在仓库根目录运行）
script_dir = fileparts(mfilename('fullpath'));
code_dir = fullfile(script_dir, '代码');

% 需要依次运行的脚本（按学习递进顺序）
runnable = {
    'pv_single_diode_model'      % 1. 光伏单二极管模型与 I-V/P-V 曲线
    'INC_algorithm'              % 2. 电导增量法
    'PandO_algorithm'            % 3. 扰动观察法（基础版）
    'perturb_and_observe_basic'  % 4. 扰动观察法（曲线实现）
    'global_scan_pv'             % 5. 全局扫描+精跟（局部遮阴多峰）
    'mppt_adaptive_slope_boost'  % 6. 变步长自适应 P&O + Boost（核心）
    'PSO_algorithm'              % 7. 粒子群算法
    'mppt_algorithm_comparison'  % 8. P&O vs PSO 对比
};

% 检测当前是否 Octave，以便用对应命令行参数
is_octave = exist('OCTAVE_VERSION', 'builtin') > 0;
if is_octave
    runner = 'octave --no-gui --eval';
else
    runner = 'matlab -batch';
end

fprintf('========== MPPT 仓库一键验证开始 ==========\n');
n_pass = 0; n_fail = 0;

for i = 1:numel(runnable)
    script = runnable{i};
    fpath = fullfile(code_dir, [script, '.m']);
    if ~exist(fpath, 'file')
        warning('跳过：找不到文件 %s', fpath);
        continue;
    end
    fprintf('\n----- [%d/%d] 运行 %s -----\n', i, numel(runnable), script);

    % 在独立子进程中运行脚本，输出重定向到临时文件，避免终端噪音
    logf = fullfile(tempdir, ['run_' script '.log']);
    cmd = sprintf('%s "run(''%s'')" > "%s" 2>&1', runner, fpath, logf);
    [status, ~] = system(cmd);

    if status == 0
        n_pass = n_pass + 1;
        fprintf('     ✓ %s 运行完成\n', script);
    else
        n_fail = n_fail + 1;
        fprintf('     ✗ %s 运行失败（exit=%d）\n', script, status);
        % 打印前几行错误便于定位
        fid = fopen(logf, 'r');
        if fid > 0
            c = 0;
            while ~feof(fid) && c < 6
                line = fgetl(fid);
                if ischar(line)
                    fprintf('        | %s\n', line);
                end
                c = c + 1;
            end
            fclose(fid);
        end
    end
end

fprintf('\n========== 验证结束：通过 %d / %d，失败 %d ==========\n', ...
        n_pass, n_pass + n_fail, n_fail);
fprintf('如需查看单个脚本的完整输出，请查看 %srun_*.log\n', tempdir);