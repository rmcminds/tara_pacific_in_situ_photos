%% requires that you follow instructions at the following page to install listdlg2, except use my custom listdlg3 script https://www.mathworks.com/matlabcentral/answers/97420-how-do-i-specify-the-location-of-a-listdlg
rehash toolboxcache
clear all

close all

yesno = [true(),false(),NaN];
standardshots = ["standard_wideangle" "standard_closeup" "very_wide" "other"];
borers = ["bivalve" "Spirobranchus" "Tridacna" "boringurchin" "other_polychaete" "sponge" "other"];
algae = ["Halimeda" "Turbinaria" "Dictyota" "Lobophora" "cca" "Galaxaura" "Sargassum" "other"];
targets = ["Porites","Pocillopora","Millepora","other"];
bleached = ["normal" "light" "bleached" "dark"];
othernotable = ["no" "other" "gallcrabs" "ascidians"];
poc_morphotypes = ["grandis" "elegans" "meandrina" "verrucosa" "damicornis" "acuta" "capitata" "other"];
por_morphotypes = ["lobata" "lutea" "columnar" "knobby" "other"];
mil_morphotypes = ["platyphylla" "dichotoma" "other"];
quadrats = ["large_thin_white","small_fat_grey","large_mid_white","hammer","other","none"];
quadrat_sizes = [21.75,19^2,20.75^2,28.5*10]; %% in cm2
quadrat_sizes_size = size(quadrat_sizes);

colnames = ["filename" "isntjunk" "hastarget" "widevsclose" "target" "morphotype" "quadrat" "parrotfish_scars" "bivalve" "Spirobranchus" "Tridacna" "boringurchin" "other_polychaete" "sponge" "otherborers" "Halimeda" "Turbinaria" "Dictyota" "Lobophora" "cca" "Galaxaura" "Sargassum" "otheralgae" "bleached" "unhealthy" "turf_over" "sediment" "colsize_quadratio" "colsize_cm2" "porites_lumps" "porites_ridges" "pigmentation" "trematodiasis" "annotationsession" "annotator" "gallcrabs" "ascidians" "notes"];
coltypes = ["string","logical","logical","string","string","string","string","logical","logical","logical","logical","logical","logical","logical","string","logical","logical","logical","logical","logical","logical","logical","string","string","logical","logical","logical","double","double","logical","logical","logical","logical","string","string","logical","logical","string"];
varsize = size(colnames);

imgpath = uigetdir(); %% define folder of input images
cd(imgpath) %% move to that folder
names = dir('*CORAL*.jpg'); %% put the names of all the jpgs in the folder into an object (careful when expanding to full dataset because there could be variable file types?)
names2 = string({names.name}'); %% extract a simple list of names from the object, transpose it, and convert from character to string
                
annotationsession = datestr(clock,'yyyymmddHHMM');
annotator = string(inputdlg('What are your initials?'));
prompt = {'Do you want to take quantitative measurements?'};
doingquantitative = yesno(listdlg('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no'}));

prompt = {'Are you adding to existing annotations?'};
if yesno(listdlg('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no'}))
    [temp,existingannotations] = uigetfile('*','select existing annotation file');
    existingannotations = [existingannotations,temp];
    opts = detectImportOptions(existingannotations, 'Delimiter', '\t');
    varNames = opts.VariableNames; % variable names
    if length(varNames) ~= length(coltypes)
        uiwait(msgbox('Old annotation fields do not match new ones. Output file will only contain new annotations', 'note', 'modal'));
        oldresults = readtable(existingannotations);
        names_rand = string.empty;
        for i = 1:length(names2)
            if not(any(strcmpi(extractBefore(names2(i),'.'), extractBefore(oldresults.filename,'.'))))
                  names_rand = [names_rand; names2(i)];
            end
        end
        names_rand = names_rand(randperm(length(names_rand))); %% randomize the order
        results = table('Size', [0 varsize(2)], 'VariableTypes', coltypes, 'VariableNames', colnames);
        resstartsize = 0;
    else
        varTypes = coltypes;
        varTypes(strcmpi('logical', varTypes)) = 'double';
        opts = setvartype(opts,varNames,varTypes);
        oldresults = readtable(existingannotations, opts);
        names_rand = string.empty;
        for i = 1:length(names2)
            if not(any(strcmpi(extractBefore(names2(i),'.'), extractBefore(oldresults.filename,'.'))))
                  names_rand = [names_rand; names2(i)];
            end
        end
        names_rand = names_rand(randperm(length(names_rand))); %% randomize the order
        results = oldresults;
        resstartsize = length(results.filename);
    end
    standardshots = unique([standardshots, rmmissing(results.widevsclose)'], 'stable');
    targets = unique([targets, rmmissing(results.target)'], 'stable');
    poc_morphotypes = unique([poc_morphotypes, rmmissing(results(strcmp(results.target, 'Pocillopora'),:).morphotype)'], 'stable');
    por_morphotypes = unique([por_morphotypes, rmmissing(results(strcmp(results.target, 'Porites'),:).morphotype)'], 'stable');
    mil_morphotypes = unique([mil_morphotypes, rmmissing(results(strcmp(results.target, 'Millepora'),:).morphotype)'], 'stable');
    quadrats = unique([quadrats, rmmissing(results.quadrat)'], 'stable'); %% there might be a fourth or even fifth quadrat to watch out for
else
    results = table('Size', [0 varsize(2)], 'VariableTypes', coltypes, 'VariableNames', colnames);
    names_rand = names2(randperm(length(names2))); %% randomize the order
    resstartsize = 0;
end % if adding to preexisting annotations

h = msgbox([' ' char(int2str(length(names_rand))) ' photos to go!']);
uiwait(h)
%% start looping through photos. use the cancel button on any prompt to either end the session or redo the active photo from the beginning (technically any error will lead you to this choice)
i = resstartsize + 1;
while i < (resstartsize + length(names_rand))
    try
        close all
                       
        results = [results; array2table(nan(1,varsize(2)), 'VariableNames', colnames)]; %% append a new row to the results table that corresponds to this photo
        results.filename(i) = names_rand(i - resstartsize);
        results.annotationsession(i) = annotationsession;
        results.annotator(i) = annotator;
                       
        currentimage = imread(char(results.filename(i))); %% open the image
        imshow(currentimage)
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'units','normalized','innerposition',[0 0 1 1]);
        

        %% 1)
        prompt = {'Is this a useful photo?'};
        results.isntjunk(i) = yesno(listdlg3('PromptString',prompt,'ListString',{'yes','no'},'position',[600 0 0 0]));
        if not(results.isntjunk(i))
            prompt = {'Anything notable?'};
            if yesno(listdlg3('PromptString',prompt,'ListString',{'yes','no'},'position',[600 0 0 0]))
                results.notes(i) = inputdlg('Describe what''s notable');
            end
            close all
            i = i + 1;
            continue
        end

        %% 2)
        prompt = {'Is there a target colony in the photo?'};
        results.hastarget(i) = yesno(listdlg3('PromptString',prompt,'ListString',{'yes','no'},'position',[600 0 0 0]));
        if not(results.hastarget(i))
            prompt = {'Anything notable?'};
            if yesno(listdlg3('PromptString',prompt,'ListString',{'yes','no'},'position',[600 0 0 0]))
                results.notes(i) = inputdlg('Describe what''s notable');
            end
            close all
            i = i + 1;
            continue
        end

        %% 2.5)
        prompt = {'Is this a \"standard\" wideangle or closeup shot?'}; % standard meaning when we were doing just one of each, with quadrat in the wideangle and showing the entire colony, and no scale but showing polyp details in the closeup
        results.widevsclose(i) = standardshots(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',standardshots,'position',[600 0 0 0]));
        if strcmpi(results.widevsclose(i), "other")
            results.widevsclose(i) = inputdlg('Describe what''s nonstandard about the photo');
            standardshots = [standardshots, results.widevsclose(i)];
        end

        %% 2.5.5)
        prompt = {'Which quadrat is in the photo?'};
        whichQuad = listdlg3('PromptString',prompt,'ListString',quadrats,'position',[600 0 0 0]);
        results.quadrat(i) = quadrats(whichQuad);
        if strcmpi(results.quadrat(i), "other")
            results.quadrat(i) = string(inputdlg('Create a quadrat identifier'));
            quadrats = [quadrats, results.quadrat(i)];
        end

        %% 5)
        prompt = {'Are there parrotfish scars on the target?'};
        results.parrotfish_scars(i) = yesno(listdlg3('PromptString',prompt,'ListString',{'yes','no','unsure/maybe'},'position',[600 0 0 0]));
                       
        prompt = {'Are there boring organisms growing in the target?'};
        results.bivalve(i) = false();
        results.Spirobranchus(i) = false();
        results.Tridacna(i) = false();
        results.boringurchin(i) = false();
        results.other_polychaete(i) = false();
        results.sponge(i) = false();
        results.otherborers(i) = "none";
        if yesno(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no'},'position',[600 0 0 0]))
            prompt = {'Which borers? (you can select multiple using command)'};
            currborers = borers(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',borers,'position',[600 0 0 0]));
            if any(strcmpi(currborers, "bivalve"))
                results.bivalve(i) = true();
            end
            if any(strcmpi(currborers, "Spirobranchus"))
                results.Spirobranchus(i) = true();
            end
            if any(strcmpi(currborers, "Tridacna"))
                results.Tridacna(i) = true();
            end
            if any(strcmpi(currborers, "urchin"))
                results.boringurchin(i) = true();
            end
            if any(strcmpi(currborers, "other_polychaete"))
                results.other_polychaete(i) = true();
            end
            if any(strcmpi(currborers, "sponge"))
                results.sponge(i) = true();
            end
            if any(strcmpi(currborers, "other"))
                results.otherborers(i) = inputdlg('What''re the other borers?');
            end
        end
                       
        prompt = {'Is the target in direct contact with algae?'};
        results.Halimeda(i) = false();
        results.Turbinaria(i) = false();
        results.Dictyota(i) = false();
        results.Lobophora(i) = false();
        results.cca(i) = false();
        results.Sargassum(i) = false();
        results.Galaxaura(i) = false();
        results.otheralgae(i) = "none";
        if yesno(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no'},'position',[600 0 0 0]))
            prompt = {'Which algae? (you can select multiple holding command)'};
            curralg = algae(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',algae,'position',[600 0 0 0]));
            if any(strcmpi(curralg, "Halimeda"))
                results.Halimeda(i) = true();
            end
            if any(strcmpi(curralg, "Turbinaria"))
                results.Turbinaria(i) = true();
            end
            if any(strcmpi(curralg, "Dictyota"))
                results.Dictyota(i) = true();
            end
            if any(strcmpi(curralg, "Lobophora"))
                results.Lobophora(i) = true();
            end
            if any(strcmpi(curralg, "cca"))
                results.cca(i) = true();
            end
            if any(strcmpi(curralg, "Sargassum"))
                results.Sargassum(i) = true();
            end
            if any(strcmpi(curralg, "Galaxaura"))
                results.Galaxaura(i) = true();
            end
            if any(strcmpi(curralg, "other"))
                results.otheralgae(i) = inputdlg('What''re the other algae?');
            end
        end
                       
        prompt = {'Is the target in contact with sediment?'};
        results.sediment(i) = yesno(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no','unsure/maybe'},'position',[600 0 0 0]));
        prompt = {'Does the target have patches of turf overgrowth?'};
        results.turf_over(i) = yesno(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no','unsure/maybe'},'position',[600 0 0 0]));
        prompt = {'Is the target bleached?'};
        results.bleached(i) = bleached(listdlg3('PromptString',prompt,'ListString',bleached,'position',[600 0 0 0]));
        prompt = {'Does the target have a pigmentation response? (usually pink)'};
        results.pigmentation(i) = yesno(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no','unsure/maybe'},'position',[600 0 0 0]));
        if not(isnan(results.pigmentation(i)))
            if results.pigmentation(i)
                prompt = {'Is it trematodiasis'};
                results.trematodiasis(i) = yesno(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no','unsure/maybe'},'position',[600 0 0 0]));
            end
        end
        prompt = {'Is the target otherwise unhealthy or recently losing tissue?'};
        results.unhealthy(i) = yesno(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',{'yes','no','unsure/maybe'},'position',[600 0 0 0]));

        %% 3)
        prompt = {'Which target genus is the colony?'};
        results.target(i) = targets(listdlg3('PromptString',prompt,'ListString',targets,'position',[600 0 0 0]));
        if strcmpi(results.target(i), "other")
            results.target(i) = inputdlg('What other target is it?');
            targets = [targets, results.target(i)];
        end

        prompt = {'Anything else notable? (you can select multiple holding command)'};
        results.gallcrabs(i) = false();
        results.ascidians(i) = false();
        results.notes(i) = 'none';
        notablenow = othernotable(listdlg3('PromptString',prompt,'ListSize',[300,300],'ListString',othernotable,'position',[600 0 0 0]));
        if not(strcmpi(notablenow, "no"))
            if any(strcmpi(notablenow, "gallcrabs"))
                results.gallcrabs(i) = true();
            end
            if any(strcmpi(notablenow, "ascidians"))
                results.ascidians(i) = true();
            end
            if any(strcmpi(notablenow, "other"))
                results.notes(i) = inputdlg('Describe what''s notable');
            end
        end

        if strcmpi(results.target(i), "Pocillopora")

            %% 4)
            prompt = {'Which morphotype is the colony?'};
            results.morphotype(i) = poc_morphotypes(listdlg3('PromptString',prompt,'ListString',poc_morphotypes, 'position', [600 0 0 0]));
            if strcmpi(results.morphotype(i), "other")
                results.morphotype(i) = string(inputdlg('What other morphotype is it?'));
                poc_morphotypes = [poc_morphotypes, results.morphotype(i)];
            end

            if doingquantitative

                %% 6)
                title('click on the inside corners of the quadrat, clockwise from the bottom left')
                movingPoints = [ginput(1); ginput(1); ginput(1); ginput(1)]; %% identify the inner corners of the quadrat
                fixedPoints = [0 1000; 1000 1000; 1000 0; 0 0]; %% not sure why these are the proper coordinates; I'm guessing MatLab must have the vertical axis with positive numbers downward?? Since the first coordinate is positive on the vertical axis instead of zero?
                tform = fitgeotrans(movingPoints, fixedPoints, 'Projective'); %% find transformation to squish the points from the original photo into the 'fixedPoints' defined above.
                RA = imref2d([1000 1000], [0 1000], [0 1000]); %% not entirely sure
                [I,r] = imwarp(currentimage, tform, 'OutputView', RA); %% transform the photo

                I = rgb2gray(I); %% go to greyscale. this discards info - would be nice to use all color channels somehow
                [I,temp] = histeq(I); %% standardize brightness and contrast. important for consistent segmentation
                imshow(I)

                [L,N] = superpixels(I,100000); %% sort of a rough segmentation; finds regions of low contrast and merges pixels to smooth things out. I actually found this an important step to sort of standardize the inputs - since each photo will have a different number of pixels inside the cropped quadrat, the pixels will be at different realworld scales, which can introduce a bias. this step normalizes each image so that the region inside the quadrat always has the same number of 'pixels'.
                outputImage = zeros(size(I),'like',I);
                idx = label2idx(L);
                numRows = size(I,1);
                numCols = size(I,2);
                for labelVal = 1:N
                    redIdx = idx{labelVal};
                    outputImage(redIdx) = mean(I(redIdx));
                end
                IS = outputImage;
                imshow(IS)

                gmag = imgradient(IS); %% next few sections from https://www.mathworks.com/help/images/marker-controlled-watershed-segmentation.html
                imshow(gmag,[])

                se = strel('disk',20);

                Ie = imerode(I,se);
                Iobr = imreconstruct(Ie,IS);
                imshow(Iobr)


                Iobrd = imdilate(Iobr,se);
                Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
                Iobrcbr = imcomplement(Iobrcbr);
                imshow(Iobrcbr)
                title('Opening-Closing by Reconstruction')

                fgm = imregionalmax(Iobrcbr);
                imshow(fgm)
                I2 = labeloverlay(I,fgm);
                imshow(I2)

                se2 = strel(ones(5,5));
                fgm2 = imclose(fgm,se2);
                fgm3 = imerode(fgm2,se2);
                fgm4 = bwareaopen(fgm3,20);
                I3 = labeloverlay(I,fgm4);
                imshow(I3)

                T = adaptthresh(Iobrcbr,0.85,'NeighborhoodSize', 2*floor(size(Iobrcbr)/32)+1, 'Statistic', 'median');
                bw = imbinarize(Iobrcbr, T); %% this is dependent on image scale! a cropped image with fewer pixels will have more sensitivity relative to pixels; the same sensitivity if the cropped image has the same realworld scale
                %imshow(bw)
                D = bwdist(bw);
                DL = watershed(D);
                bgm = DL == 0;
                %imshow(bgm)
                %title('Watershed Ridge Lines)')

                gmag2 = imimposemin(gmag, bgm | fgm4);
                L = watershed(gmag2);

                %labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm4;
                %I4 = labeloverlay(I,labels);
                %imshow(I4)
                %title('Markers and Object Boundaries Superimposed on Original Image')

                imshow(I)
                set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
                set(gcf, 'units', 'normalized', 'innerposition', [0 0 1 1]);
                hold on
                Lrgb = label2rgb(L,'jet','w','shuffle');
                himage = imshow(Lrgb);
                himage.AlphaData = 0.2;
                title('Left-click on all the segments that correspond to target coral branches. Right click to stop; ESC to quit')

                particproperties = regionprops(L);
                particname=[];
                for j = 1:size(particproperties,1)
                    particcentroids = [particproperties(j).Centroid];
                    particname = [particname, L(round(particcentroids(2)), round(particcentroids(1)))];
                    plot(round(particcentroids(1)), round(particcentroids(2)), 'ko')
                    text(round(particcentroids(1)), round(particcentroids(2)), num2str(double(particname(j)))); %% label each segmented region
                end

                totregion = zeros(1000,1000);
                moreimage = imshow(label2rgb(ones(1000,1000),'jet',[0 0 0],'shuffle'));
                moreimage.AlphaData = 0;
                datmat2 = uint8.empty;
                button = 1;
                while button==1

                    [xi,yi,button] = ginput(1);
                    xi = round(xi);
                    yi = round(yi);
                    selected = L(round(yi),round(xi));
                    region = L==selected;
                    totregion = totregion+region;
                    moreimage.AlphaData = 0.5 * totregion;

                    if(button == 3) %If right button is pressed finish digitising
                    %    h = msgbox('Stop selecting! ');
                    %    uiwait(h)
                    elseif(button == 27) %If ESC is pressed
                        h = msgbox(' Stop! Digitalization Interrupted......');
                        uiwait(h)
                        close
                        clear
                        break;
                    else % if left button is used
                        plot(xi,yi,'b+')
                        datmat2 = [datmat2; xi yi selected];

                    end %if button3

                end %while button1 %% select the relevant regions

                totregion = totregion~=0; %%binary mask; if a pixel belonged to at least one of the clicked regions it will now belong to this superregion
                
                results.colsize_quadratio(i) = sum(sum(totregion)) / (1000*1000); %% since the mask is a binary array, the sum of all elements is the number of pixels that belong to the mask. we define the number of pixels for the image above as 1000*1000, so this quantity is the percentage of pixels in the cropped photo that corresponds to the coral colony
                if whichQuad <= quadrat_sizes_size(2)
                    results.colsize_cm2(i) = results.colsize_quadratio(i) * quadrat_sizes(whichQuad); %% multiply the percent of the image corresponding to coral by the real-world scale of the image to get the real-world size of the colony. this is a container that accumulates the results for each image as the loop goes on
                end

                close all
            end % if doingquantitative
        elseif strcmpi(results.target(i), "Porites")
            %% 4)
            prompt = {'Which morphotype is the colony?'};
            results.morphotype(i) = por_morphotypes(listdlg3('PromptString',prompt,'ListString',por_morphotypes,'position', [600 0 0 0]));
            if strcmpi(results.morphotype(i), "other")
                results.morphotype(i) = string(inputdlg('What other morphotype is it?'));
                por_morphotypes = [por_morphotypes, results.morphotype(i)];
            end

            prompt = {'Does the colony have round lumpy features?'};
            results.porites_lumps(i) = yesno(listdlg3('PromptString',prompt,'ListString',{'yes','no','unsure/maybe'},'position', [600 0 0 0]));
            prompt = {'Are there linear ridges on the surface of the colony?'};
            results.porites_ridges(i) = yesno(listdlg3('PromptString',prompt,'ListString',{'yes','no','unsure/maybe'},'position', [600 0 0 0]));
        elseif strcmpi(results.target(i), "Millepora")
            %% 4)
            prompt = {'Which morphotype is the colony?'};
            results.morphotype(i) = mil_morphotypes(listdlg3('PromptString',prompt,'ListString',mil_morphotypes,'position', [600 0 0 0]));
            if strcmpi(results.morphotype(i), "other")
                results.morphotype(i) = string(inputdlg('What other morphotype is it?'));
                mil_morphotypes = [mil_morphotypes, results.morphotype(i)];
            end
        end % if strcmpi(results.target(i), "Pocillopora")

        i = i + 1;
    catch
        results(i,:) = [];
        prompt = {'Do you want to end the session or redo this photo?'};
        if yesno(listdlg('PromptString', prompt, 'ListString', {'end session', 'redo photo'}))
            break
        end
    end % try
    writetable(results, ['results_' annotationsession '.txt'], 'Delimiter', '\t')
end % while i < (resstartsize + length(names_rand))

close all

writetable(results, ['results_' annotationsession '.txt'], 'Delimiter', '\t')
                       
if i == (resstartsize + length(names_rand))
  h = msgbox(['Congratulations! You''re done! (I think)']);
  uiwait(h)
end

%% consider using this as template for creating input for transfer learning - eg. identify individual branches and use them as input for neural network; click on all branches to have the entire colony segmented, and use that as input to NN, and then use 'the rest' as input to neural network separately. the outputs from the NNs can then be treated as four separate datasets in factor analysis; microstructure, macrostructure, environment within quadrat, and environment outside of quadrat. of course also use this to get specific quantities - colony size (area of all branches normalized to quadrat[total colony pixels per subsetted image]), branch density (number of segmented regions / colony size), branch size (not sure how to semi-automatrically determine because more of a branch will be visible if density is lower, so density and length are confounded), rugosity (region circumference / region area; will depend on size and density of verrucae and on sub-branching patterns; might want to tackle those measurements more directly by re-doing segmentation of individual branches)
