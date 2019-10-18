# automated-nctools

1. function [status,msg,msgID,mk_path_fin]=rootdir(root,dir_lvl,dir_names)

This function creates directories relative to specified root directory

Inputs:
 - root: root directory where the directories are to be created 
 - dir_lvl: level of directory, relative to root directory (row vector)
 - dir_names: name of each correspoding directory (row vector)

 Outputs:
 - for all other options (k-nearest points)
     1. status -> Folder creation status indicating whether the attempt to
          create the folder is successful, returned as 0 or 1. 
          If the attempt to create the folder is successful or the folder 
          already exists,then the value of status is 1. Otherwise, the
          value is 0.
     2. msg -> Error message, returned as a character vector. If an error
          or warning occurs, msg contains the message text of the error or
          warning. Otherwise, msg is empty, ''.
	  3. msgID -> Error message identifier, returned as a character vector.
          If an error or warning occurs, msgID contains the message identifier
          of the error or warning. Otherwise, msgID is empty, ''.
	  4. mk_path_fin -> the yielded full file paths.

to do list:
 - automatic report generation
 - improve processing speed via mexing or not ->pre-allocate all cells
 and matrixes (too large paths of directory cannot be created because
 there is a limit and actually there is no point of such paths
 but I want to try to make it faster)
 - GUI
Coded by Voulanas Dimitrios on 2019-05-12
Copyright (c) 2019 Dimitrios Voulanas. 
http://https://github.com/dvoulanas
e-mail: dvoulanas@yahoo.com
All rights reserved.
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are 
met:
	* Redistributions of source code must retain the above copyright 
    notice, this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright 
   notice, this list of conditions and the following disclaimer in the 
    documentation and/or other materials provided with the distribution.
	* The names of its contributor must not be used to endorse or promote products 
    derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY Dimitrios Voulanas ''AS IS'' AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL Dimitrios Voulanas BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
