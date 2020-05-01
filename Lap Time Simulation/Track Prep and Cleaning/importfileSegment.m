function [VarName1, VarName2, VarName3] = importfileSegment(workbookFile, sheetName, dataLines)
%IMPORTFILE1 Import data from a spreadsheet
%  [VARNAME1, VARNAME2, VARNAME3] = IMPORTFILE1(FILE) reads data from
%  the first worksheet in the Microsoft Excel spreadsheet file named
%  FILE.  Returns the data as column vectors.
%
%  [VARNAME1, VARNAME2, VARNAME3] = IMPORTFILE1(FILE, SHEET) reads from
%  the specified worksheet.
%
%  [VARNAME1, VARNAME2, VARNAME3] = IMPORTFILE1(FILE, SHEET, DATALINES)
%  reads from the specified worksheet for the specified row interval(s).
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  [x, y, segment] = importfileSegment("M:\projects\Dhillon-00AEROUNCCSD1\SD 2\Lap Time Simulation\tracks\Cleaned\winstonSegmented.xlsx.xlsx", "Sheet1", [2, 952]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 09-Mar-2020 22:52:46

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 952];
    
end

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 3);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":C" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
opts.VariableTypes = ["double", "double", "double"];

% Import the data
tbl = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":C" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    tbl = [tbl; tb]; %#ok<AGROW>
end

%% Convert to output type
VarName1 = tbl.VarName1;
VarName2 = tbl.VarName2;
VarName3 = tbl.VarName3;
end