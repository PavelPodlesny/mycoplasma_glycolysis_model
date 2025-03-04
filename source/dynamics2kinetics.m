

function [kinetics] = dynamics2kinetics(nva)

        % nva stands for name-value arguments
        % Command line-example:
        % dynamics2kinetics(predictor="ATP", dynamics_data_path=data_path, init_cond_path=ic_path, output_dir=output_dir, fitType='custom');
        arguments
                nva.predictor          char % name of experimentally measured metabolite
                nva.predictor_type     char % 'product' or 'reagent'
                nva.dynamics_data_path char
                nva.output_dir         char
                nva.init_cond_path     char
                nva.fitType            char = 'rat11' %'rat11', 'cubicspline', 'custom', 'piecewise'
                nva.plot_graph         logical = false

        end

        % read experimental data
        dynamics = getDynamicData(nva);
        init_cond = getInitialCond(nva);

        % get unique IDs
        rctIDs = unique(dynamics.GroupID);
        % in the file rctIDs are duplicated because there are multiple time points in one assay (assay == unique rctID)

        % allocate memory for kinetics data
        sz = [numel(rctIDs) 8]; % numel = number of elements
        varTypes = [repmat({'string'},[1 3]), repmat({'double'},[1 5])];
        varNames = {'ID', 'InitCond','Variable','VariableValue','V0','V0_err','BioRep','TechRep'}; % [conc.] = mkM, [time] = sec --> [V0] = mkM/sec
        kinetics = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

        % for each assay (a technical repeat as a rule) calculate V0 as derivative of fitting curve at time=0
        pageIs1st = true; % for plotting diagnostical graphs
        for i=1:length(rctIDs)
                % rctID (GroupID from dynamic data table): "assay_R_B_T", where assay -- name of varied metabolite,
                % R/B/T -- number of a reaction/bio repeat/technical repeat
                ID2split = split(rctIDs(i),"_"); % "f6p_1_1_1" --> ["f6p", "1", "1", "1"]
                variable = upper(ID2split(1)); 
                
                variable_value = init_cond{strcmp(init_cond.ID,rctIDs(i)) & strcmp(init_cond.Name,variable), ["PropertyValue"]};
                ic = zipInitialCond(init_cond, rctIDs(i)); % from tabular format to one-string format ("m1=c1:m2=c2:m3=c3")
                
                techRep = dynamics(dynamics.GroupID==rctIDs(i), ["Time", nva.predictor]);

                curve = struct;
                curve.name = rctIDs(i);
                curve.x = techRep{:, ["Time"]};
                curve.y = techRep{:, [nva.predictor]};
                [curve.y_trans, v0_sign] = transConc(curve.y);
                [v0, curve.v0_err, fitObj] = estimateSlope(curve, nva.fitType, 0); % 3rd arg = (x0==0)
                % v0 is calculated in a way that it is positive always
                if (strcmp(nva.predictor_type,'reagent'))
                        if(v0_sign<0)
                                curve.v0 = v0; 
                        else %v0_sign>=0
                                curve.v0 = v0*(-1);
                        end

                elseif (strcmp(nva.predictor_type,'product'))
                        if(v0_sign<0)
                                curve.v0 = v0*(-1); 
                        else % v0_sign>=0
                                curve.v0 = v0;
                        end
                end
                %curve.v0 = v0*v0_sign; %?
                
                % add kinetics to a table, plot fitted curve
                kinetics(i,:) = {strjoin(ID2split(1:3),"_"), ic, variable, ...
                                 variable_value, curve.v0, curve.v0_err, str2num(ID2split(3)), str2num(ID2split(4))}; %!!! abs(curve.v0)?
                %plot_path = [nva.output_dir filesep 'fitting_quality.pdf'];
                if nva.plot_graph
                        plotFit2Dynamics(curve, nva.fitType, fitObj, [nva.output_dir filesep 'fitting_quality_' nva.fitType '.pdf'], pageIs1st);
                end
                
                pageIs1st = false;
        end
        kinetics = sortrows(kinetics, ["ID"]);
        writetable(kinetics, [nva.output_dir filesep 'kinetics_' nva.fitType '.csv'], ...
                        'FileType', 'text', ...
                        'Delimiter', ',');

        % aggregate technical repeats to bio-repeats
        % rctID = "STAGE_R_B.T", R -- reaction, B -- bio-rep, T -- tech-rep
        % cutTechRepID(rctID) = "STAGE_R_B"
        %subIDs = unique(cutTechRepID(kinetics.rctID, delim));   
end

function tbl = getDynamicData(nva)
        
        opts = detectImportOptions(nva.dynamics_data_path, ...
                                   'FileType', 'text', ...
                                   'TextType', 'string', ...
                                   'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        opts.VariableTypes = [{'string'}, repmat({'double'}, 1, length(opts.VariableNames) - 1)];
        tbl = readtable(nva.dynamics_data_path, opts);
end

function tbl = getInitialCond(nva)
        opts = detectImportOptions(nva.init_cond_path, ...
                                   'FileType', 'text', ...
                                   'TextType', 'string', ...
                                   'Delimiter', '\t', ...
                                   'ReadVariableNames', true, 'ReadRowNames', false);
        opts.VariableTypes = [repmat({'string'},[1 4]), {'double'}];
        tbl = readtable(nva.init_cond_path, opts);
end

function ic = zipInitialCond(init_cond, rctID)

       subset = init_cond(strcmp(init_cond.ID, rctID), ["Name", "PropertyValue"]);
       ic = strings(1, height(subset));
       for i = 1:height(subset)
               ic(i) = sprintf('%s=%d', subset.Name{i}, subset.PropertyValue(i));
       end
       ic = strjoin(ic, ':');
        
end

function [slope, err, fitObj] = estimateSlope(data, fitType, x0)

        err = NaN; %dummy values
        
        if strcmp(fitType, 'piecewise')

                lmspecs = ["constant" repmat("linear",1,length(data.x)-3); ...
                            repmat("linear",1,length(data.x)-3) "constant" ];
                fitObj = struct;
                min_rmse = Inf;
                start_brkpt = 2; %start break point
                for i=0:(length(data.x)-3)
                        brkpt = start_brkpt + i;
                        x1 = data.x(1:brkpt);
                        y1 = data.y_trans(1:brkpt);
                        lm1 = fitlm(x1, y1, lmspecs(1,i+1));
                        res1 = lm1.Residuals.Raw';
                        
                        x2 = data.x(brkpt:end);
                        y2 = data.y_trans(brkpt:end);
                        lm2 = fitlm(x2, y2, lmspecs(2,i+1));                        
                        res2 = lm2.Residuals.Raw';

                        res = [res1(1:end-1) absmax([res1(end), res2(1)]) res2(2:end)];
                        rmse = lm1.RMSE + lm2.RMSE;
                        %disp([brkpt rmse]);
                        if (rmse<min_rmse)
                                min_rmse = rmse;
                                fitObj.brkpt = brkpt;
                                fitObj.lm1 = lm1;
                                fitObj.lm2 = lm2;
                                fitObj.res = res;
                                fitObj.RMSE = rmse;
                        end
                end
                if (length(data.x)<=5) %few points available
                        lm1 = fitlm(data.x, data.y_trans);
                        %disp([0 lm1.RMSE]);
                        lm2 = fitlm(data.x(1:end-1), data.y_trans(1:end-1));
                        %disp([-1 lm2.RMSE]);
                        [rmse, i] = min([lm1.RMSE lm2.RMSE min_rmse]);
                        switch i
                                case 1
                                        fitObj.brkpt = 0;
                                        fitObj.lm1 = lm1;
                                        fitObj.lm2 = NaN;
                                        fitObj.res = lm1.Residuals.Raw';
                                case 2
                                        fitObj.brkpt = -1;
                                        fitObj.lm1 = lm2;
                                        fitObj.lm2 = NaN;
                                        fitObj.res = [lm2.Residuals.Raw' 0];
                                otherwise
                                        %% pass
                                        
                        end
                        fitObj.RMSE = rmse;
                end
                
                slope = fitObj.lm1.Coefficients.Estimate(2);
                err = fitObj.lm1.Coefficients.SE(2);
                
        else
                if strcmp(fitType, 'custom')
                        fo = fitoptions('Method','NonlinearLeastSquares',...
                                'Lower',[0, 0],... % [min(0, data.y_trans(end)), 0]
                                'Upper',[max(data.y_trans)*2, data.x(end)],... % [max(0, data.y_trans(end)), data.x(end)]
                                'StartPoint',[data.y_trans(end), data.x(end)/2]);
                        fitType = fittype('a*x/(x+b)', 'coefficients', ["a" "b"], 'options', fo);
                end
                fitObj = fit(data.x, data.y_trans, fitType);
                %slope = differentiate(fitObj, x0);
                cvals = coeffvalues(fitObj);
                cnames = coeffnames(fitObj);
                %disp(cnames); %debug!!!
                %disp(cvals);  %debug!!!
                slope = cvals(1)/cvals(2);
                
                ci = confint(fitObj);
                n_samples = 1000;
                d_samples = zeros(1, n_samples);
                
                for i = 1:n_samples
                    a = rand * (ci(2,1) - ci(1,1)) + ci(1,1);
                    b = rand * (ci(2,2) - ci(1,2)) + ci(1,2);
                    d_samples(i) = a/b;
                end
                
                err = std(d_samples);
        end
        
end

function plotFit2Dynamics(data, fittype, fitObj, plot_path, pageIs1st)
        %settings
        alpha = 0.95;
        fs = 10;
        x_offset = 50;
        y_offset = 100;
        
        fig = figure('Visible','off');
        % transformed y + PI
        subplot(2,2,1);
        ax = gca; ax.FontSize = 8;
        %plot(fitObj, "b-", data.x, data.y, "ko", "predobs", alpha); hold on;
        if strcmp(fittype, 'piecewise')
                plot(ax, fitObj.lm1); hold on;
                if (fitObj.brkpt~=-1)
                        plot(ax, fitObj.lm2); hold on;
                end
        else
                % visualize v0
                slope = data.v0;
                intercept = data.y_trans(1) - sign(slope)*slope*data.x(1);
                x = linspace(data.x(1), data.x(1)+x_offset, 100);
                y = sign(slope)*slope*x + intercept;
                % prediction interval
                pred = predint(fitObj, data.x, alpha ,'observation','off');
                plot(fitObj, "b-", data.x, data.y_trans, "ko"); hold on;
                plot(data.x, pred,'m--'); hold on;
                plot(x, y, 'r-');
        end
        xlim([min(data.x) max(data.x)+x_offset]); ylim([min(data.y_trans) max(data.y_trans)+y_offset]);
        %xlim([0 (max(data.x)+10)]); ylim([0 (min(data.y)+10)]);
        title([string(['Fit: PI/CI ' num2str(100*alpha) '%, trans. dynamics'])], ...
               'Interpreter', 'none', 'FontSize', fs);
        ylabel('Units', 'FontSize', fs); xlabel('Time, sec.', 'FontSize', fs);
        legend off;

        subplot(2,2,2);
        ax = gca; ax.FontSize = 8;
        plot(data.x, data.y, 'bo');
        title( ['original dynamics'], ...
               'Interpreter', 'none', 'FontSize', fs);
        ylabel('Conc., mkM', 'FontSize', fs); xlabel('Time, sec.', 'FontSize', fs);
        legend off;

        subplot(2,2,3);
        ax = gca; ax.FontSize = 8;
        if strcmp(fittype, 'piecewise')
                plot(data.x, fitObj.res, "bo"); hold on;
                plot(data.x, zeros(1, length(data.x)), "r-");
        else
                plot(fitObj, data.x, data.y_trans, "residuals");
        end 
        ylabel('Residuals, mkM', 'FontSize', fs); xlabel('Time, sec.', 'FontSize', fs);
        legend off;

        subplot(2,2,4)
        ax = gca; ax.FontSize = 8;
        plot(ones(1,5), 1:5, 'wo');
        if strcmp(fittype, 'piecewise')
                txt = { string(['assay ID: ' char(data.name)]), ...
                        ['fitType: ' fittype], ...
                        ['V0 = ' num2str(round(data.v0,2))], ...
                        ['V0_err = ' num2str(round(data.v0_err,2))], ...
                        ['RMSE = ' num2str(round(fitObj.RMSE,2))], ...
                        ['break pt: ' num2str(fitObj.brkpt)]};
        else
                txt = { string(['assay ID: ' char(data.name)]), ...
                        ['fitType: ' fittype], ...
                        ['V0 = ' num2str(round(data.v0,2))], ...
                        ['V0_err = ' num2str(round(data.v0_err,2))]};
        end
        
        text(0.5,3,txt, ...
            'FontSize', 12, ...
            'Interpreter', 'none');
        
        if pageIs1st
                %disp(plot_path);
                exportgraphics(fig, plot_path, 'BackgroundColor', 'none', 'ContentType', 'vector');
        else
                exportgraphics(fig, plot_path, 'Append', true, 'BackgroundColor', ...
                                'none', 'ContentType', 'vector');
        end 
        close(fig);
        
end

function [y_, reflect] = transConc(y)
        eps = 1e-10;
        shift = y(1); % y(t=0)
        reflect = sign(eps+min(y - shift));
        y_ = reflect*(y - shift);
end

function amax = absmax(v)
        [~, i] = max(abs(v));
        amax = v(i);
end

function amin = absmin(v)
        [~, i] = min(abs(v));
        amin = v(i);
end