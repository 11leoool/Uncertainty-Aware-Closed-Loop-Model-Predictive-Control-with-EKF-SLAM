% PAIRED_REANALYSIS_V2 - corrected paired analysis of the shared-seed studies.
%
% Supersedes paired_reanalysis.m, which had four defects identified in review:
%
%  (1) INVALID p-VALUES. It bootstrapped around the OBSERVED paired mean and
%      measured how often that distribution crossed zero. That inverts a
%      confidence interval; it does not generate a null distribution. p-values
%      here come from a paired SIGN-FLIP PERMUTATION test (the exact null for a
%      symmetric paired difference), while the bootstrap is used only for CIs.
%      The old floor value "5e-05" merely meant "no resample crossed zero".
%
%  (2) WRONG DIRECTION LABEL. A single "negative is better" header is false for
%      clearance and waypoint metrics, where larger is better. Each metric now
%      carries its own direction.
%
%  (3) wp_reached IS NOT WAYPOINT COMPLETION. It is the executive's stage index,
%      which also advances on timeout, so every mode returns 3. True-frame
%      waypoint attainment is recomputed here from the stored true trajectory
%      (a waypoint counts as visited if the TRUE pose passes within 0.15 m).
%
%  (4) PATH LENGTH IS NOT A CLEAN EFFICIENCY COST. Simulations continue after a
%      collision, so a comparator that collides in 44% of trials contributes
%      post-collision path. Path length is therefore reported over ALL trials
%      (labelled as simulated path) AND over the collision-free subset of pairs.
%
% Three further repairs after a second review:
%
%  (5) SEED. The resampling was previously unseeded, so the intervals and
%      permutation p-values were not reproducible run to run. rng(31,'twister')
%      is set below and is the only source of randomness in this script.
%
%  (6) COLLISION MULTIPLICITY. The four exact McNemar tests are now one declared
%      family and are Holm-corrected together, reported after both studies.
%      Previously only the continuous outcomes were corrected.
%
%  (7) COLLISION-FREE PATH IS A SELECTED SUBGROUP. Conditioning on BOTH
%      controllers avoiding a collision selects on an outcome the treatment
%      itself affects, so the surviving pairs are not exchangeable and a
%      hypothesis test on them carries no protected error rate. That row is now
%      EXPLORATORY: paired CI only, no p-value, excluded from every family.
%
% Multiplicity: within each study the confirmatory continuous outcomes form one
% Holm family (declared here, before the results are read). The four collision
% tests form a second family, pooled across both studies. The exploratory row
% belongs to neither.
%
% Run: matlab -batch "paired_reanalysis_v2"

rng(31,'twister');                 % the only randomness in this script
B = 20000; NPERM = 20000; WP_TOL = 0.15;
fid = fopen('paired_reanalysis_v2.txt','w');
pr = @(varargin) fprintf(fid,varargin{:});
pr('=== Corrected paired re-analysis (no new simulation) ===\n');
pr('Seed: rng(31,''twister'') - resampling is reproducible.\n');
pr('CIs: paired bootstrap over trials (%d). p-values: paired sign-flip\n', B);
pr('permutation (%d), which is the correct null for a paired difference.\n', NPERM);
pr('Holm correction within each study over its confirmatory continuous family,\n');
pr('and across the four collision tests as a separate declared family.\n');
pr('Direction is stated per metric.\n\n');

STUD = { 'um_obstacle_results.mat','Static sensed obstacle',{'odom_fixed','slam_fixed'},'slam_cov'; ...
         'um_dynamic_results.mat', 'Dynamic obstacle',      {'static','cv_fixed'},      'cv_cov' };

PC = []; LC = {};                  % collision family, pooled across studies

for s = 1:size(STUD,1)
    fn = STUD{s,1};
    if ~isfile(fn), pr('--- %s missing ---\n\n', fn); continue; end
    S = load(fn); R = S.results; ref = STUD{s,4}; wps = S.cfg.wps;
    pr('================ %s ================\n', STUD{s,2});
    P = []; L = {};
    for ci = 1:numel(STUD{s,3})
        cmp = STUD{s,3}{ci};
        if ~isfield(R,cmp), continue; end
        a = R.(ref); b = R.(cmp); n = numel(a.collided);
        ca = logical(a.collided(:)); cb = logical(b.collided(:));
        b10 = sum(ca & ~cb); b01 = sum(~ca & cb); nd = b10+b01;
        if nd>0, pmc = min(1,2*binocdf(min(b01,b10),nd,0.5)); else, pmc = 1; end
        rr = (sum(cb)-sum(ca))/n;                     % paired risk reduction
        bsr = zeros(B,1);
        for t=1:B, idx=randi(n,n,1); bsr(t)=(sum(cb(idx))-sum(ca(idx)))/n; end
        rlo=quantile(bsr,0.025); rhi=quantile(bsr,0.975);
        pr('\n-- %s vs %s (n=%d) --\n', ref, cmp, n);
        pr('  collisions %d/%d vs %d/%d | discordant %d:%d | exact McNemar p=%.3g\n', ...
            sum(ca),n,sum(cb),n,b10,b01,pmc);
        pr('  paired risk reduction %+.1f pp [%.1f, %.1f]\n', 100*rr, 100*rlo, 100*rhi);
        if b10==0, pr('  no trial contained a %s-only collision.\n', ref); end
        PC(end+1)=pmc; LC{end+1}=sprintf('%s: %s vs %s', STUD{s,2}, ref, cmp); %#ok<AGROW>
        % ---- true-frame waypoints ----
        wa = truewp(a.sample, wps, WP_TOL); wb = truewp(b.sample, wps, WP_TOL);
        % ---- continuous outcomes ----
        okpair = ~ca & ~cb;                            % collision-free pairs
        % cols: name | ref vals | cmp vals | direction | in Holm family | test it?
        M = { 'terminal error', a.term_true_ref(:), b.term_true_ref(:), 'lower better', true, true;
              'min clearance',  a.min_clear(:),     b.min_clear(:),     'HIGHER better', true, true;
              'simulated path (all)', a.path_len(:), b.path_len(:),     'lower better', true, true;
              'path, collision-free pairs [EXPLORATORY]', a.path_len(okpair), b.path_len(okpair), 'lower better', false, false;
              'true waypoints /3', wa(:), wb(:),                        'HIGHER better', true, true };
        for k = 1:size(M,1)
            x = M{k,2}; y = M{k,3};
            if isempty(x) || numel(x)~=numel(y) || numel(x)<3, continue; end
            keep = isfinite(x) & isfinite(y);        % never test across NaNs
            if nnz(keep) < 3
                pr('  %-42s SKIPPED (insufficient finite pairs)\n', M{k,1}); continue;
            end
            x = x(keep); y = y(keep);
            d = x - y; m = numel(d);
            bs = zeros(B,1);
            for t=1:B, idx=randi(m,m,1); bs(t)=mean(d(idx)); end
            lo=quantile(bs,0.025); hi=quantile(bs,0.975);
            if ~M{k,6}
                % selected subgroup: describe it, do not test it
                pr('  %-42s %+8.4f  [%+8.4f,%+8.4f]  no test       (%s, n=%d)\n', ...
                    M{k,1}, mean(d), lo, hi, M{k,4}, m);
                continue
            end
            % paired sign-flip permutation null
            obs = abs(mean(d)); cnt = 0;
            for t=1:NPERM
                sgn = (randi(2,m,1)*2-3);
                if abs(mean(d.*sgn)) >= obs, cnt = cnt+1; end
            end
            pperm = (cnt+1)/(NPERM+1);
            pr('  %-42s %+8.4f  [%+8.4f,%+8.4f]  perm p=%.4g   (%s, n=%d)\n', ...
                M{k,1}, mean(d), lo, hi, pperm, M{k,4}, m);
            if M{k,5}, P(end+1)=pperm; L{end+1}=sprintf('%s | %s',cmp,M{k,1}); end %#ok<AGROW>
        end
    end
    % Holm within this study's declared confirmatory family
    if ~isempty(P)
        [sp,ord]=sort(P); mm=numel(P); adj=zeros(1,mm); run=0;
        for i=1:mm, run=max(run,(mm-i+1)*sp(i)); adj(i)=min(1,run); end
        pr('\n  -- Holm within study (family of %d confirmatory continuous tests) --\n', mm);
        for i=1:mm
            pr('    %-42s perm p=%-9.4g Holm=%-9.4g %s\n', L{ord(i)}, sp(i), adj(i), ...
                tern(adj(i)<0.05,'survives','ns'));
        end
        pr('  (the exploratory collision-free path row is excluded from this family)\n');
    end
    pr('\n');
end

% ---- Holm across the four collision tests (one declared family) ----
if ~isempty(PC)
    [sp,ord]=sort(PC); mm=numel(PC); adj=zeros(1,mm); run=0;
    for i=1:mm, run=max(run,(mm-i+1)*sp(i)); adj(i)=min(1,run); end
    pr('================ Collision family (Holm across all %d tests) ================\n', mm);
    for i=1:mm
        pr('  %-52s McNemar p=%-11.4g Holm=%-11.4g %s\n', LC{ord(i)}, sp(i), adj(i), ...
            tern(adj(i)<0.05,'survives','ns'));
    end
    pr('\n');
end

% ---------------- randomized geometry: TWO-WAY bootstrap ----------------
gf = '../randomized-geometry/um_randgeom_v2_results.mat';
if isfile(gf)
    G = load(gf); modes = cellstr(string(G.modes));
    NG = double(G.NG); MS = double(G.M_EVAL);
    % rebuild [NG x MS] collision matrices per mode from per-trial records
    C = zeros(NG,MS,numel(modes));
    for i=1:numel(G.TRIAL)
        t=G.TRIAL(i); mi=find(strcmp(modes,t.mode));
        C(t.geom,t.seed,mi)=double(t.collided_cont);
    end
    ic = find(strcmp(modes,'slam_cov'));
    pr('================ Randomized geometry v2 ================\n');
    pr('The SAME 20 noise realizations are reused across all 30 geometries, so the\n');
    pr('design is CROSSED (geometry x seed), not merely clustered within geometry.\n');
    pr('Intervals below resample geometries AND seeds independently (two-way), which\n');
    pr('removes the conditioning on the fixed noise panel that a geometry-only\n');
    pr('bootstrap imposes.\n\n');
    for j=1:numel(modes)
        if j==ic, continue; end
        D = C(:,:,ic) - C(:,:,j);
        pt = mean(D(:));
        bs = zeros(B,1);
        for t=1:B
            gi = randi(NG,NG,1); si = randi(MS,MS,1);
            sub = D(gi,si); bs(t) = mean(sub(:));
        end
        lo=quantile(bs,0.025); hi=quantile(bs,0.975);
        pr('  slam_cov vs %-12s rate diff %+7.4f  two-way 95%% CI [%+7.4f,%+7.4f]\n', ...
            modes{j}, pt, lo, hi);
    end
    pr('\n  Note: a no-difference result here is absence of evidence, not evidence of\n');
    pr('  equivalence; no equivalence margin was specified in advance.\n\n');
end
fclose(fid); type paired_reanalysis_v2.txt
fprintf('\nPAIRED_V2_DONE\n');

function w = truewp(samp, wps, tol)
    % the static trials store the true trajectory as traj_true, the dynamic
    % trials as traj; accept either rather than silently returning NaN
    n = numel(samp); w = zeros(n,1);
    for i=1:n
        if isempty(samp{i}), w(i)=NaN; continue; end
        if isfield(samp{i},'traj_true'),  T = samp{i}.traj_true(1:2,:);
        elseif isfield(samp{i},'traj'),   T = samp{i}.traj(1:2,:);
        else, w(i)=NaN; continue; end
        c = 0;
        for k=1:size(wps,2)
            if min(vecnorm(T - wps(:,k))) < tol, c = c+1; end
        end
        w(i)=c;
    end
end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
