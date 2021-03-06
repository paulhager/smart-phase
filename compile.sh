#!/bin/bash

mkdir -p ./compile
javac -cp ./Phase/src/src/smartPhase ./Phase/src/src/smartPhase/*.java -d ./compile/ -classpath ./Phase/libraries/picard.jar:./Phase/libraries/commons-cli-1.4.jar 
cp ./Phase/libraries/*.jar ./compile/
cd compile/
jar xf picard.jar
jar xf commons-cli-1.4.jar
rm picard.jar
rm commons-cli-1.4.jar
jar -cvfm smartPhase.jar ../Phase/META-INF/MANIFEST.MF .
