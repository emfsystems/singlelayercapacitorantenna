classdef SingleLayerCapacitorAntenna_new
    % SingleLayerCapacitorAntenna_new
    %
    % This class models a single-layer capacitor antenna array that uses
    % microstrip-based feed lines, offers multi-port feeding, supports
    % optional compound sub-array definitions, visualizes geometry on a PCB,
    % and computes far-field patterns. It also includes methods to handle
    % reception, bandwidth/impedance analyses, noise calculations, and
    % time-domain impulse responses.

    properties
        %% --- Core Antenna Properties ---
        NumberOfElements      % Total count of capacitor-based elements in the array
        Length                % Physical length of a single capacitor element (m)
        Thickness             % Thickness of the radiator (m)
        Alength               % Array's overall length dimension (m)
        Awidth                % Array's overall width dimension (m)
        Frequency             % Main or center operating frequency (Hz)
        Voltage               % Drive/Feed voltage in transmit mode (V)
        Capacitance           % Capacitance (F) for each single-layer capacitor element
        Permeability          % Magnetic permeability (H/m) used in some field calculations
        Distance              % Far-field distance assumption (m)
        Inductance            % Lumped inductance (H) for the antenna elements
        Resistance            % Lumped parasitic resistance (Ohms) for capacitor losses
        SourceImpedance       % Source (generator) impedance (Ohms)

        %% --- General Metadata & Solver Options ---
        AntennaType           % String describing the antenna type/family
        Description           % Brief textual description of design or usage
        Application           % Possible use-case or domain (e.g., HF comm)
        MultiPort             % Boolean indicating if multiple feeds/ports are separately driven

        %% --- Frequency & Pattern Data ---
        FrequencyVector       % Vector of frequencies for sampling (Hz)
        FrequencyResponse     % Relative frequency response values (same size as FrequencyVector)
        AntennaElement        % `phased.CustomAntennaElement` storing pattern data

        % Raw pattern data for the antenna (indexed by elevation, azimuth, frequency)
        MagnitudeDataRaw      % Normalized magnitude pattern data: size [nEl x nAz x nFreq]
        PhaseDataRaw          % Phase data (radians): size [nEl x nAz x nFreq]

        %% --- Matching Network (optional, can be empty) ---
        MatchingType          % Type of matching approach used (e.g., 'none', 'lumped', etc.)
        MatchingComponents    % Struct listing matching components if needed

        %% --- Current & Impedance Tracking ---
        ElementCurrents       % Stores the final solved currents for all elements (Tx mode)

        %% --- Array Configuration ---
        ConfigurationType     % 'single', 'rectangular', or 'compound'
        NumRows               % Number of rows (for rectangular arrays)
        NumCols               % Number of columns (for rectangular arrays)
        CompoundArrayDefinitions   % Struct array specifying sub-arrays if in 'compound' config

        %% --- Geometry & Feed Arrangements ---
        CapacitorCoords       % Nx3 array of (x,y,z) coordinates for each capacitor element
        FeedPoints            % Struct array describing each feed point (coords, voltage, microstrip info, etc.)
        AzimuthAngles         % Vector of azimuth angles (deg) used for pattern calculations
        ElevationAngles       % Vector of elevation angles (deg) used for pattern calculations

        %% --- Visualization Toggles ---
        DrawFeedLines   = true;
        ColorSubArrays  = true;
        LabelSubArrays  = true;
        CenterCompound  = true;  % If true, the entire compound array is shifted so its bounding box is centered

        %% --- Transmission-Line Phase Offsets ---
        LinePhaseOffsets      % Nx1 vector of additional phase offsets (rad) from feed line effects

        %% --- Per-element Feed Line Lengths ---
        ElementFeedLineLengths   % Nx1 array giving feed line lengths from feed point to each element

        %% --- Reception/Imaging Properties ---
        OperatingMode = 'transmit';      % 'transmit', 'receive', or 'sensor'
        ReceiverImpedance = 50;          % Typical input impedance of the receiver front-end
        LNA_Gain = 0;                    % Low-noise amplifier gain in dB
        LNA_NoiseFigure = 1;             % LNA noise figure in dB
        ReceiverNoiseTemperature = 290;  % Receiver noise temperature in Kelvin

        ReceivedVoltages            % Stores feed voltages for the array when in receive mode
        GroundPlaneEnabled  = true; % If true, simulate ground-plane reflection
        GroundPlaneZ        = 0;    % Z-plane coordinate for ground plane
        ReflectionPhaseFlip = false;% If true, flip phase upon ground-plane reflection
    end

    methods
        %% ------------------------- CONSTRUCTOR -------------------------
        function obj = SingleLayerCapacitorAntenna_new(...
                numberOfElements, lengthVal, thicknessVal, AlengthVal, AwidthVal, ...
                frequencyVal, voltageVal, capacitanceVal, permeabilityVal, distanceVal, ...
                inductanceVal, resistanceVal, sourceImpedanceVal, varargin)
            % SingleLayerCapacitorAntenna_new  Constructor method
            %
            %   obj = SingleLayerCapacitorAntenna_new(N, len, thick, A_len, A_wid,
            %       freq, volt, cap, mu, dist, L, R, Z0, ... 'Name', Value, ...)
            %
            %   This constructor configures fundamental properties of the
            %   single-layer capacitor antenna array. The geometry and feed
            %   arrangement are then defined, followed by a precomputation
            %   of radiation patterns.

            % Check for required inputs
            if nargin < 12
                error('Insufficient number of required inputs to constructor.');
            end
            % Default to 50 Ohms source if not provided
            if isempty(sourceImpedanceVal)
                sourceImpedanceVal = 50;
            end

            % Basic validation
            validateattributes(numberOfElements, {'numeric'}, {'positive','integer'});
            obj.NumberOfElements = numberOfElements;
            obj.Length           = lengthVal;
            obj.Thickness        = thicknessVal;
            obj.Alength          = AlengthVal;
            obj.Awidth           = AwidthVal;
            obj.Frequency        = frequencyVal;
            obj.Voltage          = voltageVal;
            obj.Capacitance      = capacitanceVal;
            obj.Permeability     = permeabilityVal;
            obj.Distance         = distanceVal;
            obj.Inductance       = inductanceVal;
            obj.Resistance       = resistanceVal;
            obj.SourceImpedance  = sourceImpedanceVal;

            % Basic metadata for clarity
            obj.AntennaType   = 'Enhanced Single-Layer Capacitor Antenna w/ Microstrip Feeds';
            obj.Description   = 'Includes microstrip feed lines, substrate info, multi-port support.';
            obj.Application   = 'Suited for HF to microwave range with line phase shifts.';
            obj.MultiPort     = false;  % By default, single-port drive

            % Default frequency sampling around center frequency
            numFreqSamples = 5;
            obj.FrequencyVector   = linspace(obj.Frequency*0.8, obj.Frequency*1.2, numFreqSamples);
            obj.FrequencyResponse = ones(size(obj.FrequencyVector));

            % Default angle coverage for patterns
            obj.AzimuthAngles   = -180:180;
            obj.ElevationAngles = -90:90;

            % Default matching approach is none
            obj.MatchingType       = 'none';
            obj.MatchingComponents = struct();

            % Parse optional name-value pairs for configuration
            p = inputParser;
            p.CaseSensitive = false;
            addParameter(p, 'ConfigurationType', 'rectangular', @(x) ischar(x) || isstring(x));
            addParameter(p, 'NumRows', 1, @(x) isnumeric(x) && x>0);
            addParameter(p, 'NumCols', numberOfElements, @(x) isnumeric(x) && x>0);
            addParameter(p, 'CompoundArrayDefinitions', [], @isstruct);
            addParameter(p, 'MultiPort', false, @islogical);
            addParameter(p, 'DrawFeedLines', true, @islogical);
            addParameter(p, 'ColorSubArrays', true, @islogical);
            addParameter(p, 'LabelSubArrays', true, @islogical);
            addParameter(p, 'CenterCompound', true, @islogical);
            parse(p, varargin{:});
            res = p.Results;

            % Assign parsed parameters to object
            obj.ConfigurationType        = lower(res.ConfigurationType);
            obj.NumRows                  = res.NumRows;
            obj.NumCols                  = res.NumCols;
            obj.CompoundArrayDefinitions = res.CompoundArrayDefinitions;
            obj.MultiPort                = res.MultiPort;
            obj.DrawFeedLines            = res.DrawFeedLines;
            obj.ColorSubArrays           = res.ColorSubArrays;
            obj.LabelSubArrays           = res.LabelSubArrays;
            obj.CenterCompound           = res.CenterCompound;

            % Build geometry & feed arrangement based on config
            obj = obj.defineGeometryAndFeeds();

            % Initialize line-based phase offsets
            obj.LinePhaseOffsets = zeros(obj.NumberOfElements, 1);

            % Precompute radiation pattern data for each frequency in FrequencyVector
            [magnitudePattern, phasePattern] = obj.computePattern(obj.AzimuthAngles, obj.ElevationAngles);
            obj.MagnitudeDataRaw = magnitudePattern;
            obj.PhaseDataRaw     = phasePattern;

            % Create a corresponding phased.CustomAntennaElement
            obj.AntennaElement = phased.CustomAntennaElement(...
                'FrequencyVector',   obj.FrequencyVector, ...
                'FrequencyResponse', obj.FrequencyResponse, ...
                'AzimuthAngles',     obj.AzimuthAngles, ...
                'ElevationAngles',   obj.ElevationAngles, ...
                'MagnitudePattern',  magnitudePattern, ...
                'PhasePattern',      phasePattern);
        end

        %% ------------- DEFINING GEOMETRY & FEED POINTS -------------
        function obj = defineGeometryAndFeeds(obj)
            % defineGeometryAndFeeds  Sets up the element coordinates and feed
            % configurations based on 'single', 'rectangular', or 'compound' types.

            switch obj.ConfigurationType
                case 'single'
                    % Single-element geometry is placed at origin with a single feed
                    obj.CapacitorCoords = [0, 0, 0];

                    sF.Coord   = [0, 0, 0];
                    sF.Label   = 'SingleFeed';
                    sF.Color   = 'r';
                    sF.Voltage = obj.Voltage;
                    % Example microstrip parameters
                    sF.Microstrip = struct(...
                        'Width',                0.001, ...
                        'Thickness',            0.000035, ...
                        'SubstrateHeight',      0.0016, ...
                        'EpsilonR',             4.4, ...
                        'LossTangent',          0.02, ...
                        'ConductorConductivity',5.8e7 ...
                    );
                    sF.BoardSize = [0, 0];

                    obj.FeedPoints = sF;

                case 'rectangular'
                    % Rectangular array of elements, with one feed for all
                    if obj.NumRows * obj.NumCols ~= obj.NumberOfElements
                        error('Rectangular config mismatch: NumRows*NumCols must equal NumberOfElements.');
                    end
                    obj.CapacitorCoords = obj.assignRectangularArray(...
                        obj.NumRows, obj.NumCols, obj.Alength, obj.Awidth);

                    % Choose a feed coordinate at array edge or center
                    halfWid = obj.Awidth / 2;
                    sF.Coord   = [0, halfWid, 0];
                    sF.Label   = 'RectArrayFeed';
                    sF.Color   = 'r';
                    sF.Voltage = obj.Voltage;
                    % Example microstrip params
                    sF.Microstrip = struct(...
                        'Width',                0.003, ...
                        'Thickness',            0.00005, ...
                        'SubstrateHeight',      0.001, ...
                        'EpsilonR',             4.4, ...
                        'LossTangent',          0.02, ...
                        'ConductorConductivity',5.8e7 ...
                    );
                    sF.BoardSize = [0, 0];
                    obj.FeedPoints = sF;

                case 'compound'
                    % Multiple sub-arrays, each with its own feed (common in phased arrays)
                    def = obj.CompoundArrayDefinitions;
                    if isempty(def)
                        error('Must provide CompoundArrayDefinitions for compound config.');
                    end

                    coordsAll = [];
                    feedArr   = [];
                    totalCount = 0;

                    for iSub = 1:numel(def)
                        nR   = def(iSub).NumRows;
                        nC   = def(iSub).NumCols;
                        nEl  = nR * nC;
                        totalCount = totalCount + nEl;

                        subLen = def(iSub).Alength;
                        subWid = def(iSub).Awidth;

                        % Generate local sub-array coordinates
                        subCoords = obj.assignRectangularArray(nR, nC, subLen, subWid);

                        % Apply any array-level offset specified
                        if isfield(def(iSub), 'ArrayOffset')
                            offset = def(iSub).ArrayOffset;
                        else
                            offset = [0, 0, 0];
                        end
                        subGlobal = subCoords + offset;
                        coordsAll = [coordsAll; subGlobal]; %#ok<AGROW>

                        % Build feed struct for this sub-array
                        halfWid   = subWid / 2;
                        offsetf   = [0, halfWid, 0];
                        offset    = offset + offsetf;
                        newFeed.Coord = offset;

                        if isfield(def(iSub), 'FeedOffset')
                            newFeed.Coord = newFeed.Coord + def(iSub).FeedOffset;
                            newFeed.FeedOffset = def(iSub).FeedOffset;
                        end

                        newFeed.Label = sprintf('Sub#%d', iSub);
                        if isfield(def(iSub), 'Label')
                            newFeed.Label = def(iSub).Label;
                        end

                        % Assign color
                        colCycle = {'b','g','m','c','y','k','r',[0.5, 0.5, 0],[0.8, 0.2, 0.6]};
                        newFeed.Color = colCycle{mod(iSub-1, numel(colCycle)) + 1};
                        if isfield(def(iSub), 'Color')
                            newFeed.Color = def(iSub).Color;
                        end

                        % Determine feed voltage if multi-port or single
                        if obj.MultiPort && isfield(def(iSub), 'Voltage')
                            newFeed.Voltage = def(iSub).Voltage;
                        elseif obj.MultiPort
                            newFeed.Voltage = 1;
                        else
                            newFeed.Voltage = obj.Voltage;
                        end

                        % Microstrip properties if provided
                        if isfield(def(iSub), 'Microstrip')
                            newFeed.Microstrip = def(iSub).Microstrip;
                        else
                            newFeed.Microstrip = struct(...
                                'Width',                0.003, ...
                                'Thickness',            0.00005, ...
                                'SubstrateHeight',      0.001, ...
                                'EpsilonR',             4.4, ...
                                'LossTangent',          0.02, ...
                                'ConductorConductivity',5.8e7 ...
                                );
                        end

                        % Optional board geometry
                        if isfield(def(iSub), 'BoardSize')
                            newFeed.BoardSize = def(iSub).BoardSize;
                        else
                            newFeed.BoardSize = [0, 0];
                        end

                        feedArr = [feedArr, newFeed]; %#ok<AGROW>
                    end

                    % Confirm total elements match NumberOfElements
                    if totalCount ~= obj.NumberOfElements
                        error('Compound config mismatch: sum of sub-arrays must match NumberOfElements.');
                    end

                    % Center entire arrangement if requested
                    if obj.CenterCompound
                        [coordsAll, feedArr] = obj.centerCompoundArrangement(coordsAll, feedArr);
                    end

                    % Final assignment
                    obj.CapacitorCoords = coordsAll;
                    obj.FeedPoints      = feedArr;

                otherwise
                    error('Unknown ConfigurationType: %s', obj.ConfigurationType);
            end

            % After geometry is set, compute feed-line lengths from each feed to each element
            N = obj.NumberOfElements;
            obj.ElementFeedLineLengths = zeros(N, 1);

            for elemIdx = 1:N
                [subArrIdx, ~] = obj.mapElementToSubArray(elemIdx);
                if subArrIdx == 0
                    subArrIdx = 1; % fallback if single/rect config
                end
                elemCoord = obj.CapacitorCoords(elemIdx, :);
                feedCoord = obj.FeedPoints(subArrIdx).Coord;

                dx = elemCoord(1) - feedCoord(1);
                dy = elemCoord(2) - feedCoord(2);
                dz = elemCoord(3) - feedCoord(3);

                dist = sqrt(dx^2 + dy^2 + dz^2);
                obj.ElementFeedLineLengths(elemIdx) = dist;
            end
        end

        function [coordsCentered, feedArrCentered] = centerCompoundArrangement(obj, coordsAll, feedArr)
            % centerCompoundArrangement  Translates geometry so that the bounding
            % box of all elements is centered at (0,0) in the XY-plane.

            if isempty(coordsAll)
                coordsCentered = coordsAll;
                feedArrCentered = feedArr;
                return;
            end

            minX = min(coordsAll(:,1));
            maxX = max(coordsAll(:,1));
            minY = min(coordsAll(:,2));
            maxY = max(coordsAll(:,2));
            midX = (minX + maxX) / 2;
            midY = (minY + maxY) / 2;

            coordsCentered = coordsAll;
            coordsCentered(:,1) = coordsAll(:,1) - midX;
            coordsCentered(:,2) = coordsAll(:,2) - midY;

            feedArrCentered = feedArr;
            for k = 1:numel(feedArr)
                feedArrCentered(k).Coord(1) = feedArr(k).Coord(1) - midX;
                feedArrCentered(k).Coord(2) = feedArr(k).Coord(2) - midY;
            end
        end

        function coords = assignRectangularArray(obj, nRows, nCols, subLen, subWid)
            % assignRectangularArray  Generates a rectangular grid of (x,y,0)
            % coordinates representing nRows x nCols capacitors.

            % Distances between consecutive elements in each dimension
            elemSepl = subLen/(nCols - 1);
            elemSepw = subWid/(nRows - 1);
            halfLen  = subLen / 2;
            halfWid  = subWid / 2;

            [cols, rows] = meshgrid(0:(nCols-1), 0:(nRows-1));
            xPos = cols * elemSepl - halfLen;
            yPos = rows * elemSepw - halfWid;

            coords = [xPos(:), yPos(:), zeros(numel(xPos), 1)];
        end

        %% ==================== IMPEDANCE CALCULATIONS ====================
        function [Z, obj] = lineBasedImpedance(obj, freq, elemIndex)
            % lineBasedImpedance  Computes the total input impedance seen
            % by the feed for a particular element, combining:
            %   1) The "bare" single-layer capacitor element R-L-C model
            %   2) The microstrip feed line segment (if non-negligible length)

            % First compute the local R-L-C impedance
            Z_antenna = obj.antennaImpedanceSingleLayer(freq, elemIndex);

            % Identify which sub-array feed this element belongs to
            [subArrayIndex, ~] = obj.mapElementToSubArray(elemIndex);
            if subArrayIndex == 0
                subArrayIndex = 1;  % fallback
            end
            msParam = obj.FeedPoints(subArrayIndex).Microstrip;

            % Look up feed line length for this element
            feedLineLen = obj.ElementFeedLineLengths(elemIndex);

            % If feed line is extremely short, do not apply extra transformations
            if feedLineLen <= 1e-12
                Z = Z_antenna;
                return;
            end

            % Otherwise, compute the microstrip effect via rfckt classes
            [Z_in, phShift] = obj.computeMicrostripInputImpedance(Z_antenna, freq, msParam, feedLineLen);
            obj.LinePhaseOffsets(elemIndex) = phShift;  % track phase offset for reference
            Z = Z_in;
        end

        function [subIndex, localIdx] = mapElementToSubArray(obj, elemIdx)
            % mapElementToSubArray  Returns the sub-array index that the
            % given element belongs to. For single/rect, always returns (1, elemIdx).

            subIndex = 0;
            localIdx = 0;

            if ~strcmpi(obj.ConfigurationType, 'compound')
                % Non-compound => single feed
                subIndex = 1;
                localIdx = elemIdx;
                return;
            end

            % For compound arrays, check each sub-array definition
            def = obj.CompoundArrayDefinitions;
            elCounter = 1;
            for s = 1:numel(def)
                nR = def(s).NumRows;
                nC = def(s).NumCols;
                nEl = nR * nC;
                if elemIdx >= elCounter && elemIdx < elCounter + nEl
                    subIndex = s;
                    localIdx = elemIdx - elCounter + 1;
                    return;
                end
                elCounter = elCounter + nEl;
            end
        end

        function Z = antennaImpedanceSingleLayer(obj, freq, ~)
            % antennaImpedanceSingleLayer  R-L-C model for the single-layer capacitor.
            %   Ignores mutual coupling from other elements.

            omega = 2*pi*freq;
            C = obj.Capacitance;
            L = obj.Inductance;
            R = obj.Resistance;

            Z = R + 1j*(omega*L - 1/(omega*C));
        end

        function antennaLoadCkt = createAntennaLoadCircuit(R, L, C)
            % createAntennaLoadCircuit  Builds an rfckt.series circuit with R, L, and C.
            Rckt = rfckt.resistor('R', R);
            Lckt = rfckt.inductor('L', L);
            Cckt = rfckt.capacitor('C', C);
            antennaLoadCkt = rfckt.series('Ckts', {Rckt, Lckt, Cckt});
        end

        function [Zin, phaseShift] = computeMicrostripInputImpedance(obj, ~, freq, msParam, feedLineLen)
            % computeMicrostripInputImpedance  Cascades a microstrip line
            % segment with the single-layer capacitor's R-L-C circuit to find
            % the net input impedance. The function also extracts an approximate
            % phase shift from the S21 term.

            % Retrieve R, L, and C from object properties
            Rval = obj.Resistance;
            Lval = obj.Inductance;
            Cval = obj.Capacitance;

            % Create a microstrip feed line of length feedLineLen
            msLine = rfckt.microstrip( ...
                'LineLength',         feedLineLen, ...
                'Width',              msParam.Width, ...
                'Height',             msParam.SubstrateHeight, ...
                'EpsilonR',           msParam.EpsilonR, ...
                'LossTangent',        msParam.LossTangent, ...
                'Thickness',          msParam.Thickness, ...
                'SigmaCond',          msParam.ConductorConductivity ...
            );

            % Create a series R-L-C to represent the element's lumped model
            antennaLoadCkt = createAntennaLoadCircuit(Rval, Lval, Cval);

            % Cascade the microstrip line with the R-L-C load
            fullCkt = rfckt.cascade('Ckts', {msLine, antennaLoadCkt});

            % Analyze the combined 2-port circuit at the specified frequency
            freqGHz = freq;
            analyze(fullCkt, freqGHz);

            % Extract the S-parameters at the single frequency point
            S = fullCkt.AnalyzedResult.S_Parameters;  % 2x2x1
            S11 = S(1,1,1);
            S21 = S(2,1,1);

            % Usually 50 ohms reference, but can vary
            zref = fullCkt.AnalyzedResult.Z0;
            if numel(zref) > 1
                zref = zref(1);
            end

            % Convert S11 to input impedance
            Zin = zref * (1 + S11)/(1 - S11);

            % Approximate net phase shift from angle of S21
            phaseShift = angle(S21);
        end

        function Zmn = mutualImpedance(obj, freq, m, n)
            % mutualImpedance  Placeholder approximation for mutual coupling
            % between elements m and n. For a more rigorous approach,
            % incorporate full electromagnetic models.

            c = physconst('LightSpeed');
            lambda = c / freq;

            d = norm(obj.CapacitorCoords(m, :) - obj.CapacitorCoords(n, :));
            if d == 0
                Zmn = obj.lineBasedImpedance(freq, m);
            else
                % A simplistic free-space approach, ignoring real conduction/dielectric losses
                Zmn = (obj.Permeability * exp(-1j * 2*pi * d / lambda)) / (2*pi*d);
            end
        end

        %% ================== PATTERN & CURRENT SOLUTIONS =================
        function [magnitudePattern, phasePattern] = computePattern(obj, az, el)
            % computePattern  Evaluates a custom far-field pattern by summing
            % contributions from each element’s current (including optional
            % ground-plane reflection). Normalizes to a max of 1 in power pattern.

            azRad = deg2rad(az);
            elRad = deg2rad(el);
            nAz   = numel(az);
            nEl   = numel(el);
            nFreq = numel(obj.FrequencyVector);

            magnitudePattern = zeros(nEl, nAz, nFreq);
            phasePattern     = zeros(nEl, nAz, nFreq);

            c       = physconst('LightSpeed');
            mu      = obj.Permeability;
            R       = obj.Distance;               % far-field spherical radius
            epinot  = 8.854187817e-12;            % permittivity of free space

            [cosEl, sinEl] = deal(cos(elRad(:)), sin(elRad(:)));
            zgrid = R * sinEl * ones(1, nAz);     % approximate observation points

            for fIdx = 1:nFreq
                freq = obj.FrequencyVector(fIdx);
                k    = 2*pi*freq/c;
                om   = 2*pi*freq;

                % Solve element currents at this frequency
                I_elem = obj.currentDistribution(freq);

                % Initialize accumulators for H-like (MF) and E fields
                patternData = zeros(nEl, nAz);
                EData       = zeros(nEl, nAz);

                % Sum contributions from each element
                for idx = 1:obj.NumberOfElements
                    I_e = I_elem(idx);
                    I_m = abs(I_e);
                    I_Phase = angle(I_e);

                    % Coordinates of the element
                    xc = obj.CapacitorCoords(idx, 1);
                    yc = obj.CapacitorCoords(idx, 2);
                    zc = obj.CapacitorCoords(idx, 3);

                    % Evaluate distance from element to far-field point
                    x = R * cosEl * cos(azRad);
                    y = R * cosEl * sin(azRad);
                    rReal = sqrt((x - xc).^2 + (y - yc).^2 + (zgrid - zc).^2);

                    phaseShiftReal = k * rReal;
                    cosAReal = (y - yc) ./ rReal;
                    aReal    = acos(cosAReal);
                    sinAReal = sin(aReal);

                    % Elementary far-field approximations
                    MF_Real = -(mu * I_m)./(2*pi*rReal) .* sinAReal .* cos(I_Phase + phaseShiftReal);
                    E_Real  =  (I_m)./(2*epinot*om*pi*rReal) .* sinAReal .* sin(I_Phase + phaseShiftReal);

                    % Optional ground-plane image element
                    if obj.GroundPlaneEnabled
                        zc_image = 2*obj.GroundPlaneZ - zc;
                        rImage = sqrt((x - xc).^2 + (y - yc).^2 + (zgrid - zc_image).^2);
                        phaseShiftImage = k * rImage;

                        cosAImg = (y - yc) ./ rImage;
                        aImg    = acos(cosAImg);
                        sinAImg = sin(aImg);

                        % Reflection magnitude (1.0 => full reflection here)
                        reflMag = 1.0;
                        if obj.ReflectionPhaseFlip
                            reflSign = -1;
                        else
                            reflSign = +1;
                        end

                        MF_Image = reflMag * (-(mu * I_m)./(2*pi*rImage)) .* sinAImg ...
                            .* cos(I_Phase + phaseShiftImage) * reflSign;
                        E_Image  = reflMag * ((I_m)./(2*epinot*om*pi*rImage)) .* sinAImg ...
                            .* sin(I_Phase + phaseShiftImage) * reflSign;

                        patternData(:, :) = patternData(:, :) + MF_Real + MF_Image;
                        EData(:, :)       = EData(:, :)       + E_Real + E_Image;
                    else
                        patternData(:, :) = patternData(:, :) + MF_Real;
                        EData(:, :)       = EData(:, :)       + E_Real;
                    end
                end

                % Zero out fields below the ground plane if enabled
                if obj.GroundPlaneEnabled
                    mask = (zgrid < obj.GroundPlaneZ);
                    patternData(mask) = 0;
                    EData(mask)       = 0;
                end

                % Combine H-like (patternData) & EData as an approximate power pattern
                powerPattern = abs(patternData.^2 + EData.^2);
                phasePattern(:,:,fIdx) = angle(patternData);

                % Normalize so maximum is 1
                maxVal = max(powerPattern(:));
                if maxVal > 0
                    powerPattern = powerPattern / maxVal;
                end
                magnitudePattern(:,:,fIdx) = powerPattern;
            end
        end

        function I_elem = currentDistribution(obj, freq)
            % currentDistribution  Solves for each element's current when driven
            % at a specified frequency in transmit mode. Considers mutual coupling
            % via a simplified mutual impedance matrix (Zmut).

            N = obj.NumberOfElements;
            Zmut = zeros(N, N);

            % Reset line-phase offsets each time if needed
            obj.LinePhaseOffsets = zeros(N, 1);

            % Build the NxN matrix of self & mutual impedances
            for m = 1:N
                for n = 1:N
                    if m == n
                        [Zmut(m, n), obj] = obj.lineBasedImpedance(freq, n);
                    else
                        Zmut(m, n) = obj.mutualImpedance(freq, m, n);
                    end
                end
            end

            % Construct the voltage vector (depends on multi-port usage)
            Vvec = zeros(N, 1);
            if ~obj.MultiPort
                % Single feed => same voltage at all elements
                Vvec(:) = obj.Voltage;
            else
                % Each sub-array feed can have distinct voltage
                def = obj.CompoundArrayDefinitions;
                if isempty(def)
                    Vvec(:) = obj.Voltage;
                else
                    elCounter = 1;
                    for s = 1:numel(def)
                        nR   = def(s).NumRows;
                        nC   = def(s).NumCols;
                        nEl  = nR * nC;
                        feedElIndex = elCounter;
                        for r = 1:nEl
                            Vvec(feedElIndex) = obj.FeedPoints(s).Voltage;
                            feedElIndex = feedElIndex + 1;
                        end
                        elCounter = elCounter + nEl;
                    end
                end
            end

            % Solve Zmut * I_elem = Vvec
            I_elem = Zmut \ Vvec;
            obj.ElementCurrents = I_elem;  % store for reference
        end

        function [mag2D, phase2D] = safePatternCall(obj, freq)
            % safePatternCall  Calls built-in pattern(...) on the internal
            % phased.CustomAntennaElement. If the result collapses to 1D,
            % revert to stored raw data.

            az = obj.AzimuthAngles;
            el = obj.ElevationAngles;

            [~, mag, ph] = pattern(obj.AntennaElement, freq, az, el, 'Type','efield');
            s = size(mag);
            if s(1) == 1 || s(2) == 1
                % Dimension collapse => fallback on precomputed data
                %warning('Dimension collapse in pattern() => using raw data fallback...');
                mag2D = obj.getRawPatternSlice(freq);

                [~, idx] = min(abs(obj.FrequencyVector - freq));
                if ~isempty(obj.PhaseDataRaw)
                    phase2D = obj.PhaseDataRaw(:,:,idx);
                else
                    phase2D = zeros(size(mag2D));
                end
            else
                mag2D   = mag;
                phase2D = ph;
            end
        end

        %% ==================== RECEPTION METHODS ====================
        function feedVoltages = computeReceivedVoltages(obj, freq, incField)
            % computeReceivedVoltages  Computes the voltages at each feed
            % in receive mode, given an incoming plane wave direction &
            % amplitude.

            if nargin < 2, freq = obj.Frequency; end
            if ~strcmpi(obj.OperatingMode,'receive') && ~strcmpi(obj.OperatingMode,'sensor')
                warning('computeReceivedVoltages called but OperatingMode is not receive/sensor.');
            end

            numFeeds = numel(obj.FeedPoints);
            feedVoltages = zeros(numFeeds,1);

            % Retrieve the receive pattern magnitude for the specified incField angles
            pData = obj.computeReceivePattern(freq);
            azInd = find(obj.AzimuthAngles == incField.Az,1);
            elInd = find(obj.ElevationAngles == incField.El,1);
            if isempty(azInd) || isempty(elInd)
                % If exact angle not found, pick the closest
                [~, azInd] = min(abs(obj.AzimuthAngles - incField.Az));
                [~, elInd] = min(abs(obj.ElevationAngles - incField.El));
            end
            elementFactor = pData.magnitude(elInd, azInd);

            % Convert LNA gain from dB to a linear voltage multiplier
            linearGain = 10^(obj.LNA_Gain/20);

            % Map each element to its sub-array feed index
            subArrayIdxForElement = zeros(obj.NumberOfElements,1);
            if strcmpi(obj.ConfigurationType,'compound')
                def = obj.CompoundArrayDefinitions;
                elCounter = 1;
                for s = 1:numel(def)
                    nR = def(s).NumRows;
                    nC = def(s).NumCols;
                    nEl = nR * nC;
                    for e = elCounter:(elCounter + nEl - 1)
                        subArrayIdxForElement(e) = s;
                    end
                    elCounter = elCounter + nEl;
                end
            else
                subArrayIdxForElement(:) = 1;
            end

            % Compute induced voltages per element, apply loading by the receiver
            for elemIdx = 1:obj.NumberOfElements
                subArr = subArrayIdxForElement(elemIdx);
                % Approx open-circuit element voltage from E-field
                V_oc_elem = incField.Amplitude * elementFactor;

                Zant = obj.lineBasedImpedance(freq, elemIdx);
                Zrx  = obj.ReceiverImpedance;
                % Simple voltage-divider for loaded voltage
                Z_tot  = (Zant * Zrx) / (Zant + Zrx);
                ratio  = Z_tot / Zant;
                V_loaded = V_oc_elem * ratio;

                % Accumulate into the feed's total
                feedVoltages(subArr) = feedVoltages(subArr) + V_loaded;
            end

            % Finally, apply LNA gain to each feed's voltage
            feedVoltages = feedVoltages * linearGain;
            obj.ReceivedVoltages = feedVoltages;
        end

        function patternData = computeReceivePattern(obj, freq)
            % computeReceivePattern  Gets or derives the pattern used
            % for receive calculations at a specific frequency.
            if nargin < 2
                freq = obj.Frequency;
            end

            [mag2D, phase2D] = obj.safePatternCall(freq);
            patternData.magnitude = mag2D;
            patternData.phase     = phase2D;
        end

        function outSignal = simulateReception(obj, incField, freq)
            % simulateReception  Applies computeReceivedVoltages and
            % returns the feed voltages for a given plane wave incField.
            if nargin < 3
                freq = obj.Frequency;
            end
            feedV = obj.computeReceivedVoltages(freq, incField);
            outSignal = feedV;
        end

        function rawData = performImagingScan(obj, azAngles, elAngles, freqRange, fieldAmplitude)
            % performImagingScan  Sweeps over given az/el angles and frequency
            % to produce feed voltages, primarily for imaging or direction scanning.

            if ~strcmpi(obj.OperatingMode, 'sensor')
                warning('performImagingScan called but OperatingMode ~= sensor.');
            end

            if nargin < 2 || isempty(azAngles)
                azAngles = -90:10:90;
            end
            if nargin < 3 || isempty(elAngles)
                elAngles = -30:10:30;
            end
            if nargin < 4 || isempty(freqRange)
                freqRange = obj.Frequency;
            end
            if nargin < 5 || isempty(fieldAmplitude)
                fieldAmplitude = 1;  % default normalized field amplitude
            end

            rawData = struct();
            idx = 1;

            % Triple-nested loop: freq -> El -> Az
            for f = freqRange
                for el = elAngles
                    for az = azAngles
                        incField.Az        = az;
                        incField.El        = el;
                        incField.Amplitude = fieldAmplitude;

                        feedV = obj.computeReceivedVoltages(f, incField);

                        rawData(idx).Frequency = f;
                        rawData(idx).Az        = az;
                        rawData(idx).El        = el;
                        rawData(idx).FeedVolt  = feedV;
                        idx = idx + 1;
                    end
                end
            end
        end

        %% ======================= VISUALIZATION =======================
        function show(obj)
            % show  Provides a 2D overhead plot of element placements
            % and feed lines, including optional board outlines.

            figure('Name', 'Antenna Geometry (Overhead w/ Boards)');
            hold on; grid on; axis equal;
            xlabel('X (m)'); ylabel('Y (m)');
            title(sprintf('Geometry: %s', obj.ConfigurationType));

            coords = obj.CapacitorCoords;
            feedArr = obj.FeedPoints;

            switch obj.ConfigurationType
                case 'single'
                    % Potentially draw a board rectangle if specified
                    if isfield(feedArr, 'BoardSize') && any(feedArr.BoardSize > 0)
                        bd = feedArr.BoardSize;
                        rectangle('Position',...
                            [-bd(1)/2, -bd(2)/2, bd(1), bd(2)],...
                            'FaceColor', [0.9, 1, 0.9],...
                            'EdgeColor', 'k',...
                            'LineWidth', 1.5);
                    end
                    plot(coords(1,1), coords(1,2), 'bo', 'MarkerSize', 6, 'LineWidth', 1.5);
                    plot(feedArr.Coord(1), feedArr.Coord(2), 'x',...
                        'Color', feedArr.Color, 'MarkerSize', 8, 'LineWidth', 2);

                    if obj.DrawFeedLines
                        line([feedArr.Coord(1), coords(1,1)], ...
                             [feedArr.Coord(2), coords(1,2)], 'Color', 'k', 'LineWidth', 1.2);
                    end

                case 'rectangular'
                    if isfield(feedArr, 'BoardSize') && any(feedArr.BoardSize > 0)
                        bd = feedArr.BoardSize;
                        rectangle('Position', ...
                            [-bd(1)/2, -bd(2)/2, bd(1), bd(2)], ...
                            'FaceColor', [0.9,1,0.9], ...
                            'EdgeColor', 'k', ...
                            'LineWidth',1.5);
                    end
                    plot(coords(:,1), coords(:,2), 'bo', 'MarkerSize', 6, 'LineWidth',1.5);
                    plot(feedArr.Coord(1), feedArr.Coord(2), 'x',...
                        'Color', feedArr.Color, 'MarkerSize',8,'LineWidth',2);

                    if obj.DrawFeedLines
                        for i = 1:size(coords,1)
                            line([feedArr.Coord(1), coords(i,1)], ...
                                 [feedArr.Coord(2), coords(i,2)], 'Color','k','LineWidth',1.0);
                        end
                    end

                case 'compound'
                    def = obj.CompoundArrayDefinitions;
                    elCounter = 1;
                    for s = 1:numel(def)
                        nR = def(s).NumRows;
                        nC = def(s).NumCols;
                        nEl = nR * nC;
                        subCoords = coords(elCounter : elCounter + nEl - 1, :);
                        elCounter = elCounter + nEl;

                        subColor = feedArr(s).Color;

                        % If board outline is given
                        if isfield(feedArr(s), 'BoardSize') && any(feedArr(s).BoardSize > 0)
                            bd = feedArr(s).BoardSize;
                            cx = feedArr(s).Coord(1);
                            cy = feedArr(s).Coord(2);
                            rectangle('Position', [cx - bd(1)/2, cy - bd(2)/2, bd(1), bd(2)], ...
                                'FaceColor',[0.9,1,0.9],'EdgeColor','k','LineWidth',1.5);
                        end
                        plot(subCoords(:,1), subCoords(:,2), 'o',...
                            'MarkerEdgeColor', subColor, 'MarkerSize',5,'LineWidth',1.0);

                        minX = min(subCoords(:,1)); maxX = max(subCoords(:,1));
                        minY = min(subCoords(:,2)); maxY = max(subCoords(:,2));
                        if obj.ColorSubArrays
                            rectangle('Position', [minX, minY, (maxX-minX), (maxY-minY)],...
                                'EdgeColor', subColor, 'LineStyle','--','LineWidth',1.5);
                        else
                            rectangle('Position', [minX, minY, (maxX-minX), (maxY-minY)],...
                                'EdgeColor','k','LineStyle','--','LineWidth',1.0);
                        end
                        if obj.LabelSubArrays
                            midX = (minX + maxX)/2;
                            midY = (minY + maxY)/2;
                            text(midX, midY, feedArr(s).Label, 'Color', subColor,...
                                'FontSize',9, 'HorizontalAlignment','center');
                        end
                        feedPt = feedArr(s).Coord;
                        plot(feedPt(1), feedPt(2), 'x', 'Color', subColor, ...
                            'MarkerSize',8,'LineWidth',2);

                        if obj.DrawFeedLines
                            for i = 1:size(subCoords,1)
                                line([feedPt(1), subCoords(i,1)],...
                                     [feedPt(2), subCoords(i,2)], 'Color','k','LineWidth',1.0);
                            end
                        end
                    end

                otherwise
                    plot(coords(:,1), coords(:,2), 'bo', 'LineWidth',1.5);
            end

            hold off;
        end

        function visualizeCurrentDistributionInteractive(obj)
            % visualizeCurrentDistributionInteractive  Demonstrates how each
            % element's phase evolves in real-time. Useful for educational
            % displays of traveling waves or phasing.

            freq = obj.Frequency;
            I_elem = obj.currentDistribution(freq);
            baseAmp   = abs(I_elem);
            basePhase = angle(I_elem) + obj.LinePhaseOffsets;

            coords = obj.CapacitorCoords;
            N      = obj.NumberOfElements;

            % Use a scaled "demo frequency" for visualization
            demoFrequency = freq * 1e-6;
            freqRad       = 2*pi*demoFrequency;

            fig = figure('Name','Interactive Transmission-Line Wave',...
                'Position',[100,100,1200,500]);

            % Left subplot: Overhead 2D scatter
            ax2D = subplot(1,2,1,'Parent',fig);
            scatter2D = scatter(ax2D, coords(:,1), coords(:,2), ...
                100*(baseAmp/max(baseAmp) + 0.1), basePhase, 'filled');
            colormap(ax2D, 'hsv');
            cb2D = colorbar(ax2D);
            cb2D.Label.String = 'Phase (rad)';
            caxis(ax2D, [-pi, pi]);
            title(ax2D, 'Overhead 2D');
            xlabel(ax2D, 'X (m)'); ylabel(ax2D, 'Y (m)');
            axis(ax2D, 'equal'); grid(ax2D, 'on');

            % Right subplot: 3D layout
            ax3D = subplot(1,2,2,'Parent',fig);
            scatter3D = scatter3(ax3D, coords(:,1), coords(:,2), coords(:,3), ...
                100*(baseAmp/max(baseAmp) + 0.1), basePhase, 'filled');
            colormap(ax3D, 'hsv');
            cb3 = colorbar(ax3D);
            cb3.Label.String = 'Phase (rad)';
            caxis(ax3D, [-pi, pi]);
            title(ax3D, '3D Perspective');
            xlabel(ax3D, 'X (m)'); ylabel(ax3D, 'Y (m)'); zlabel(ax3D, 'Z (m)');
            axis(ax3D, 'equal'); grid(ax3D, 'on');

            % Auto-scale axes
            minX = min(coords(:,1)); maxX = max(coords(:,1));
            minY = min(coords(:,2)); maxY = max(coords(:,2));
            minZ = min(coords(:,3)); maxZ = max(coords(:,3));
            marginX = 0.1 * (maxX - minX);
            marginY = 0.1 * (maxY - minY);
            marginZ = 0.1 * (maxZ - minZ);
            if (maxX <= minX), marginX = 0.1; end
            if (maxY <= minY), marginY = 0.1; end
            if (maxZ <= minZ), marginZ = 0.1; end

            xlim(ax3D, [minX - marginX, maxX + marginX]);
            ylim(ax3D, [minY - marginY, maxY + marginY]);
            zlim(ax3D, [minZ - marginZ, maxZ + marginZ]);
            view(ax3D, [-30, 20]);

            % Toggle button to animate
            btn = uicontrol(fig, 'Style','togglebutton', 'String','Play',...
                'Value',0, 'Units','Normalized', 'Position',[0.45,0.02,0.1,0.05], ...
                'BackgroundColor',[0.8,0.9,0.8], 'FontSize',10);

            % Timer-based animation
            wavePeriod = 1 / demoFrequency;
            nFramesPerPeriod = 3600;
            dt = wavePeriod / nFramesPerPeriod;
            t  = 0;

            animTimer = timer('ExecutionMode','fixedSpacing', ...
                'Period',dt, 'TimerFcn', @updateAnimation);

            function updateAnimation(~,~)
                if btn.Value == 1
                    t = t + dt;
                    rawPhase = basePhase + freqRad*t;
                    newPhase = mod(rawPhase + pi, 2*pi) - pi;
                    set(scatter2D, 'CData', newPhase);
                    set(scatter3D, 'CData', newPhase);
                    drawnow limitrate
                end
            end
            function togglePlay(src,~)
                if src.Value == 1
                    src.String = 'Pause';
                    start(animTimer);
                else
                    src.String = 'Play';
                    stop(animTimer);
                end
            end
            btn.Callback = @togglePlay;

            % Clean up timer on figure close
            fig.CloseRequestFcn = @(~,~) onClose();
            function onClose()
                stop(animTimer);
                delete(animTimer);
                delete(fig);
            end
        end

        function visualizeRadiationPattern(obj)
            % visualizeRadiationPattern  Plots the pattern in polar coords
            % at the current object frequency, sweeping all az and fixing El=0.

            freq = obj.Frequency;
            az   = obj.AzimuthAngles;
            el   = 0;

            figure;
            pattern(obj.AntennaElement, freq, az, el, ...
                'CoordinateSystem','polar','Type','powerdb');
            title('Radiation Pattern (Az sweep, El=0)');
        end

        function patternAzimuth(obj, elevationAngle)
            % patternAzimuth  Plots polar pattern vs azimuth for a specified elevation.

            freq = obj.Frequency;
            az   = obj.AzimuthAngles;
            el   = elevationAngle;

            figure;
            pattern(obj.AntennaElement, freq, az, el, ...
                'CoordinateSystem','polar','Type','powerdb');
            title(['Azimuth Pattern at El=', num2str(elevationAngle), '°']);
        end

        function patternElevation(obj, azimuthAngle)
            % patternElevation  Plots polar pattern vs elevation for a specified azimuth.

            freq = obj.Frequency;
            az   = azimuthAngle;
            el   = obj.ElevationAngles;

            figure;
            pattern(obj.AntennaElement, freq, az, el, ...
                'CoordinateSystem','polar','Type','powerdb');
            title(['Elevation Pattern at Az=', num2str(azimuthAngle), '°']);
        end

        function visualize3DPattern(obj)
            % visualize3DPattern  Displays the 3D directivity/gain pattern
            % in polar coordinates at the object’s main frequency.

            freq = obj.Frequency;
            figure;
            pattern(obj.AntennaElement, freq, 'Type','powerdb','CoordinateSystem','polar');
            title('3D Radiation Pattern');
        end

        %% ================= BASIC METRICS & INPUT/OUTPUT =================
        function D = calculateDirectivity(obj)
            % calculateDirectivity  Estimates the directivity in dB by
            % integrating the power pattern. This uses built-in pattern
            % data if available; otherwise uses stored raw fallback.

            freq = obj.Frequency;
            az   = obj.AzimuthAngles;
            el   = obj.ElevationAngles;

            [~, mag] = pattern(obj.AntennaElement, freq, az, el, 'Type','power');
            if size(mag,1) == 1
                % Fallback if dimension collapse
                %warning('Dimension collapse in calculateDirectivity => using raw fallback...');
                raw2D = obj.getRawPatternSlice(freq);
                magLinear = raw2D;
            else
                magLinear = db2pow(mag);
            end

            % Numerical integration over sphere
            deltaAz = deg2rad(az(2)-az(1));
            deltaEl = deg2rad(el(2)-el(1));
            [~, ThetaGrid] = meshgrid(deg2rad(az), deg2rad(el+90));
            Prad = sum(sum(magLinear .* sin(ThetaGrid))) * deltaAz * deltaEl;
            Umax = max(magLinear(:));
            D_lin = 4*pi * Umax / Prad;
            D = pow2db(D_lin);
        end

        function G = calculateGain(obj)
            % calculateGain  Gains the directivity + 10*log10(radiationEfficiency).

            D = obj.calculateDirectivity();
            eff = obj.calculateRadiationEfficiency();
            G = D + 10*log10(eff);
        end

        % ------------------- RADIATION RESISTANCE -------------------
        function Rr_total = calculateRadiationResistance(obj, freq)
            % calculateRadiationResistance  Computes an approximate total
            % radiation resistance for the array, summing across all elements
            % via a fraction of the capacitive reactance.

            if nargin < 2 || isempty(freq)
                freq = obj.Frequency;
            end
            omega = 2*pi*freq;
            Xc = 1/(omega * obj.Capacitance);

            % alpha is an empirical fraction for small capacitive radiators
            alpha = 1;  % Placeholder; refine from simulation or measurement
            Rrad_single_element = alpha * abs(Xc);

            % Summation logic (placeholder in code). Typically you'd do
            % Rr_total = Rrad_single_element * NumberOfElements or a more
            % advanced approach. The user-coded function call is shown:
            Rr_total = Rrad_single_element * obj.NumberOfElements;

            %Rr = Rr_total; %#ok<NASGU> 
        end

        function Rr_total = totalRadiationResistance(obj, freq)
            % totalRadiationResistance  Alternatively sums the single-element
            % Rr over the entire array (basic approximation).
            Rr_element = obj.calculateRadiationResistance(freq);
            Rr_total   = Rr_element * obj.NumberOfElements;
        end

        function Rloss = calculateLossResistance(obj, freq)
            % calculateLossResistance  Computes total loss from:
            %   (1) Per-element parasitic ESR
            %   (2) Real part of feed line impedances at each feed port

            if nargin < 2 || isempty(freq)
                freq = obj.Frequency;
            end
            Rcap = obj.Resistance;  % ESR for each element

            numFeeds = numel(obj.FeedPoints);
            totalFeedLosses = zeros(numFeeds, 1);

            subArrayIdxForElement = zeros(obj.NumberOfElements,1);

            % Identify sub-array membership for each element
            if strcmpi(obj.ConfigurationType, 'compound')
                def = obj.CompoundArrayDefinitions;
                elCounter = 1;
                for s = 1:numel(def)
                    nR = def(s).NumRows;
                    nC = def(s).NumCols;
                    nEl = nR * nC;
                    subArrayIdxForElement(elCounter : elCounter + nEl - 1) = s;
                    elCounter = elCounter + nEl;
                end
            else
                subArrayIdxForElement(:) = 1;
            end

            % For each feed, parallel combine impedances of that sub-array
            for feedIdx = 1:numFeeds
                elemIndices = find(subArrayIdxForElement == feedIdx);
                Z_feed_array = zeros(numel(elemIndices),1);

                for idx = 1:numel(elemIndices)
                    Z_feed_array(idx) = obj.lineBasedImpedance(freq, elemIndices(idx));
                end

                % Parallel combination
                Z_total_feed = 1 / sum(1./Z_feed_array);
                totalFeedLosses(feedIdx) = real(Z_total_feed);
            end

            % Sum up total conduction/dielectric losses
            totalCapLoss = Rcap * obj.NumberOfElements;
            totalFeedLoss = sum(totalFeedLosses);

            Rloss = totalCapLoss + totalFeedLoss;
        end

        function eff = calculateRadiationEfficiency(obj, freq)
            % calculateRadiationEfficiency  Aggregates each feed's total
            % radiation & loss resistances and averages if multiple feeds exist.

            if nargin < 2 || isempty(freq)
                freq = obj.Frequency;
            end
            numFeeds = numel(obj.FeedPoints);
            efficiencies = zeros(numFeeds, 1);

            subArrayIdxForElement = zeros(obj.NumberOfElements,1);
            if strcmpi(obj.ConfigurationType, 'compound')
                def = obj.CompoundArrayDefinitions;
                elCounter = 1;
                for s = 1:numFeeds
                    nR = def(s).NumRows;
                    nC = def(s).NumCols;
                    nEl = nR * nC;
                    subArrayIdxForElement(elCounter : elCounter + nEl - 1) = s;
                    elCounter = elCounter + nEl;
                end
            else
                subArrayIdxForElement(:) = 1;
            end

            for feedIdx = 1:numFeeds
                elemsInFeed = find(subArrayIdxForElement == feedIdx);
                numEl = numel(elemsInFeed);

                % Basic approach: Rr_feed = (SingleElementRr * number_of_elements_in_feed)
                elemsInFeed = find(subArrayIdxForElement == feedIdx);
                Rr_feed = obj.calculateRadiationResistance(freq) * numel(elemsInFeed);
                % --- parallel combination of impedances ---------------
                Zinv_sum = 0;                         % running sum of 1/Z for each element
                for elem = elemsInFeed'
                    Z_elem   = obj.lineBasedImpedance(freq, elem);
                    Zinv_sum = Zinv_sum + 1./Z_elem;  % accumulate reciprocal
                end
                Z_total_feed = 1 ./ Zinv_sum;         % final parallel impedance
                Rcap = obj.Resistance;
                Rloss_feed = real(Z_total_feed) + Rcap * numEl;

                efficiencies(feedIdx) = Rr_feed / (Rr_feed + Rloss_feed);
            end

            % If multiple feeds exist, take mean as overall efficiency
            eff = mean(efficiencies);
            eff = max(min(eff,1),0); % clamp in [0,1]
        end

        function estimateBandwidth(obj, freqRange, acceptableVSWR)
            % estimateBandwidth  Scans a frequency range, checking VSWR at each
            % point. The difference between min and max freq that remain below
            % the acceptable VSWR is returned as the "bandwidth".

            if nargin < 2 || isempty(freqRange)
                freqRange = linspace(obj.Frequency * 0.8, obj.Frequency * 1.2, 100);
            end
            if nargin < 3 || isempty(acceptableVSWR)
                acceptableVSWR = 2;
            end

            Z0 = obj.SourceImpedance;
            VSWR = zeros(size(freqRange));

            for idx = 1:length(freqRange)
                freq = freqRange(idx);
                Z_feed_elements = arrayfun(@(e) obj.lineBasedImpedance(freq, e), 1:obj.NumberOfElements);
                Z_total_feed = abs(1 / sum(1./Z_feed_elements));

                Gamma = (Z_total_feed - Z0) / (Z_total_feed + Z0);
                VSWR(idx) = (1 + abs(Gamma)) / (1 - abs(Gamma));
            end

            withinBW = VSWR <= acceptableVSWR;
            if ~any(withinBW)
                disp('No frequency in the specified range satisfies VSWR <= acceptable limit.');
                return;
            end

            freqBW = freqRange(withinBW);
            bandwidth = max(freqBW) - min(freqBW);
            disp(['Estimated Bandwidth: ', num2str(bandwidth/1e6), ' MHz']);
        end

        function sObj = sparameters(obj, freqRange)
            % sparameters  Builds an S-parameter object (1-port) by computing
            % reflection S11 across freqRange.

            if nargin < 2 || isempty(freqRange)
                freqRange = obj.FrequencyVector;
            end
            Z0 = obj.SourceImpedance;
            S11 = zeros(length(freqRange), 1);

            for idx = 1:length(freqRange)
                freq = freqRange(idx);
                Z_feed_elements = arrayfun(@(e) obj.lineBasedImpedance(freq, e), 1:obj.NumberOfElements);
                Z_total_feed = abs(1 / sum(1./Z_feed_elements));

                Gamma = (Z_total_feed - Z0) / (Z_total_feed + Z0);
                S11(idx) = Gamma;
            end

            sObj = obj.createSParameterObject(S11, freqRange);
        end

        function VSWR = calculateVSWR(obj, freq)
            % calculateVSWR  Computes the VSWR at a given frequency using
            % parallel-combined feed impedance and the source impedance Z0.

            if nargin < 2 || isempty(freq)
                freq = obj.Frequency;
            end
            Z0 = obj.SourceImpedance;
            Z_feed_elements = arrayfun(@(e) obj.lineBasedImpedance(freq, e), 1:obj.NumberOfElements);
            Z_total_feed = abs(1 / sum(1./Z_feed_elements));

            Gamma = (Z_total_feed - Z0)/(Z_total_feed + Z0);
            VSWR = (1 + Gamma)/(1 - Gamma);
        end

        function RL = calculateReturnLoss(obj, freq)
            % calculateReturnLoss  Computes the return loss (dB) from
            % reflection coefficient at the feed.

            if nargin < 2 || isempty(freq)
                freq = obj.Frequency;
            end
            Z0 = obj.SourceImpedance;
            Z_feed_elements = arrayfun(@(e) obj.lineBasedImpedance(freq, e), 1:obj.NumberOfElements);
            Z_total_feed = abs(1 / sum(1./Z_feed_elements));

            Gamma = (Z_total_feed - Z0)/(Z_total_feed + Z0);
            RL = -20*log10(Gamma);
        end

        function analyzePatternVsFrequency(obj)
            % analyzePatternVsFrequency  Loops through obj.FrequencyVector,
            % plotting the polar pattern at each frequency for user inspection.

            frequencies = obj.FrequencyVector;
            az = obj.AzimuthAngles;
            el = obj.ElevationAngles;
            for idx = 1:length(frequencies)
                freq = frequencies(idx);
                figure;
                pattern(obj.AntennaElement, freq, az, el, ...
                    'Type','powerdb','CoordinateSystem','polar');
                title(['Radiation Pattern at ', num2str(freq/1e9),' GHz']);
            end
        end

        function [polPatternH, polPatternV] = calculatePolarizationPattern(obj)
            % calculatePolarizationPattern  Placeholder that retrieves the
            % total power pattern and lumps it into H and V. V is zero here.

            freq = obj.Frequency;
            az   = obj.AzimuthAngles;
            el   = obj.ElevationAngles;

            [~, mag] = pattern(obj.AntennaElement, freq, az, el, 'Type','power');
            if (size(mag,1) == 1) && (length(el)>1)
                %disp('**Dimension collapse in calculatePolarizationPattern!**');
                raw2D    = obj.getRawPatternSlice(freq);
                magLinear = raw2D;
            else
                magLinear = db2pow(mag);
            end
            polPatternH = magLinear;           % hypothetical horizontal pol
            polPatternV = zeros(size(magLinear)); % placeholder vertical pol
        end

        function polEfficiency = calculatePolarizationEfficiency(obj)
            % calculatePolarizationEfficiency  Example that sums
            % "horizontal" vs total to get a fractional measure.

            [polH, polV] = obj.calculatePolarizationPattern();
            totalPower = polH + polV;
            desiredPolPower = polH;
            polEfficiency = sum(desiredPolPower(:)) / sum(totalPower(:));
        end

        function NF = calculateNoiseFigure(~)
            % calculateNoiseFigure  Placeholder for system-level noise figure (dB).
            NF = 0;
        end

        function Pn = calculateThermalNoise(~, bandwidth)
            % calculateThermalNoise  Uses k*T*B to compute noise power (W).
            if nargin < 2
                bandwidth = 1e9; % default 1 GHz
            end
            k = physconst('Boltzmann');
            T = 290;  % default temperature (K)
            Pn = k * T * bandwidth;
        end

        function exportSParameters(obj, filename, freqRange)
            % exportSParameters  Writes the single-port S-parameter data to file
            % using 'rfwrite'. The user can specify a frequency range or use
            % obj.FrequencyVector.

            if nargin < 3
                freqRange = obj.FrequencyVector;
            end
            sParams = obj.sparameters(freqRange);
            rfwrite(sParams, filename);
        end

        function exportPatternData(obj, filename, freq)
            % exportPatternData  Saves raw pattern data (magnitude & phase)
            % to a .mat file for external usage.

            if nargin < 3
                freq = obj.Frequency;
            end
            az = obj.AzimuthAngles;
            el = obj.ElevationAngles;
            [~, mag, phase] = pattern(obj.AntennaElement, freq, az, el, 'Type','efield');
            if (size(mag,1)==1) && (length(el)>1)
                %disp('**Dimension collapse in exportPatternData => fallback**');
                raw2D = obj.getRawPatternSlice(freq);
                mag   = raw2D;
                phase = zeros(size(raw2D));
            end
            save(filename, 'freq','az','el','mag','phase');
        end

        function [t, h] = impulseResponse(obj, bandwidth, numPoints)
            if nargin < 2
                bandwidth = obj.Frequency * 0.4;
            end
            if nargin < 3
                numPoints = 1024;
            end

            % Absolute frequency range (centered at operating frequency)
            freqRange_abs = linspace(obj.Frequency - bandwidth/2, obj.Frequency + bandwidth/2, numPoints);

            % Get S11 data from antenna object
            sParams = obj.sparameters(freqRange_abs);
            S11 = rfparam(sParams, 1, 1);

            % Frequency step
            deltaF = freqRange_abs(2) - freqRange_abs(1);

            % Total Sampling frequency in Hz (equal to total bandwidth)
            Fs = deltaF * numPoints;

            % Time vector symmetric around zero (suitable for TDR/impulse response)
            t = (-numPoints/2 : numPoints/2 - 1) / Fs;

            % Convert frequency-domain data (centered) to time-domain impulse response
            h = ifft(fftshift(S11), numPoints);
            h = fftshift(h); % Shift zero-frequency component to center of spectrum
        end


        function magSlice = getRawPatternSlice(obj, freqDesired)
            % getRawPatternSlice  Retrieves the magnitude pattern slice from
            % MagnitudeDataRaw at the closest frequency to freqDesired.

            if isempty(obj.MagnitudeDataRaw)
                error('MagnitudeDataRaw is empty. No stored pattern data.');
            end
            [~, idx] = min(abs(obj.FrequencyVector - freqDesired));
            magSlice = obj.MagnitudeDataRaw(:,:,idx);
        end

        function sObj = createSParameterObject(~, S11, freqRange)
            % createSParameterObject  Assembles a 1-port sparameters object
            % from a given reflection coefficient vector.

            numPorts = 1;
            sParams = zeros(numPorts,numPorts,length(freqRange));
            sParams(1,1,:) = S11;
            sObj = sparameters(sParams, freqRange);
        end
    end
end


function antennaLoadCkt = createAntennaLoadCircuit(R, L, C)
% createAntennaLoadCircuit  Builds a 2-port series RLC circuit using
% rfckt.seriesrlc from the RF Toolbox. This is used to model the R-L-C
% behavior of the single-layer capacitor element.

    antennaLoadCkt = rfckt.seriesrlc( ...
        'R', R, ...
        'L', L, ...
        'C', C ...
    );
end
