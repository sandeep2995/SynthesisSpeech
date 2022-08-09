# SynthesisSpeech
Generate cepstral coefficients from spoken speech

-------------------------------------------
Execution Instructions:
-------------------------------------------

Go to directory "Synthesis".

Set the values as desired in "config.h" file such as file name to store the speech, path of Recording Module, duration in seconds to allow recording the speech.

Inputs--->contains all inputs used for the testing 
Outputs--->contains all outputs generated for the inputs. 

Now go back to Synthesis directory (i.e. parent directory of this project)

open the program solution(Synthesis.sln) in visual studio, execute it by pressing F5, 
then it displays "Start Recording.........", you can start recording your signal. when the duration(as mentioned in "config.h" file) is over(say 4 seconds) then it will display "Stop Recording."

Now press enter to see the results.(same results will be available in out.txt file)

"cep.txt" contains the cepstral coefficient values separated by space. Each row indicates coefficents values for a frame of vowel speech

Now go to the directory Synthesis--->Synthesis.

There you can find your spoken speech in .wav as well as .txt format with the names as mentioned in "config.h" file.

scaled input speech samples can be found in "ScaledInpSpeech.txt"
