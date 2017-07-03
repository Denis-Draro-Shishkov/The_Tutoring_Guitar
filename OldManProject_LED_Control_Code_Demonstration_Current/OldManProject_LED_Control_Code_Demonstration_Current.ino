
int Multiplex1 [6] [5] = { // White, Blue, Yellow, Brown, Purple 
  {4, 3, 2, 6, 5}, // Index of Pins controlling the LEDs on the first string
  {10, 9, 8, 12, 11}, // Index of Pins controlling the LEDs on the second string
  {46, 45, 44, 48, 47}, // Index of Pins controlling the LEDs on the third string
  {40, 39, 38, 42, 41}, // Index of Pins controlling the LEDs on the fourth string
  {30, 29, 28, 32, 31}, // Index of Pins controlling the LEDs on the fifth string
  {24, 23, 22, 26, 25} // Index of Pins controlling the LEDs on the sixth string

};
int ledPin = 1;                 // LED connected to digital pin 13

int OutEnable [6] = {7, 13, 49, 43, 33, 27}; // One for each multiplexer, might want to change this later on to differentiate between outEnable and pin on. Not sure why that'd be necessary though.
int OpenLEDs [6] = {27, 35, 41, 47, 53, 7};

int i;
int Led_Number_Loop;
int StringNum_Loop;
int Time_Index_Loop;

int TestingArray[46][3] = {
  {1,0,500},
  {201,0,842},
  {3,1,500},
  {203,1,702},
  {3,2,269},
  {3,4,231},
  {203,2,269},
  {203,4,359},
  {1,5,500},
  {201,5,500},
};

void setup()
{
  Serial.begin(9600);
 
  for (int thisPin_1 = 0; thisPin_1 < 6; thisPin_1++)
  {
    for (int thisPin_2 = 0; thisPin_2 < 5; thisPin_2++)
    {
      pinMode(Multiplex1[thisPin_1] [thisPin_2], OUTPUT);
    }
  }
  for (int thisPin_3 = 0; thisPin_3 < 6; thisPin_3++)
  {
    pinMode(OutEnable[thisPin_3], OUTPUT);
  }
  for (int thisPin_4 = 0; thisPin_4 < 6; thisPin_4++)
  {
    pinMode(OpenLEDs[thisPin_4], OUTPUT);
  }
  int numRows = sizeof(TestingArray) / sizeof(TestingArray[0]);
  //Serial.print("Number of rows = ");Serial.println(numRows);

}

void loop()
{
  int numRows = sizeof(TestingArray) / sizeof(TestingArray[0]);
  Serial.print("Number of rows = "); Serial.println(numRows);
  digitalWrite(OutEnable[0], LOW); // Make sure Multiplexer is turned off before switching states.
  digitalWrite(OutEnable[1], LOW); // Make sure Multiplexer is turned off before switching states.
  digitalWrite(OutEnable[2], LOW); // Make sure Multiplexer is turned off before switching states.
  digitalWrite(OutEnable[3], LOW); // Make sure Multiplexer is turned off before switching states.
  digitalWrite(OutEnable[4], LOW); // Make sure Multiplexer is turned off before switching states.
  digitalWrite(OutEnable[5], LOW); // Make sure Multiplexer is turned off before switching states.

  for (int i = 0; i < numRows; i++)
  {
    // Serial.println(i);
    Led_Number_Loop = TestingArray[i][0];
    StringNum_Loop = TestingArray[i][1];
    Time_Index_Loop = TestingArray[i][2];
    // Serial.println(Led_Number_Loop);

    if (Led_Number_Loop < 200)
    {
      Multiplexer_Control(Led_Number_Loop, StringNum_Loop, Time_Index_Loop, OutEnable);
    }
    else
    {
      Led_Number_Loop = Led_Number_Loop - 200;
      LED_Flash(Led_Number_Loop, StringNum_Loop, Time_Index_Loop, OutEnable);
    }
  }
  digitalWrite(OutEnable[StringNum_Loop], LOW); // Turn on LED!
  while (1) {}
}

void Multiplexer_Control(int Led_Number, int StringNum, int Time_Index, int* OutEnable)
{
  if (Led_Number == 32)
  {
    digitalWrite(OpenLEDs[StringNum], 1);
    delay(Time_Index);
  }
  else
  {
    int r0 = bitRead(Led_Number, 0);   // Reads the first (least significant bit) of the LED number
    int r1 = bitRead(Led_Number, 1);   // Reads the second bit of the LED number
    int r2 = bitRead(Led_Number, 2);   // Reads the third bit of the LED number
    int r3 = bitRead(Led_Number, 3);   // Reads the fourth bit of the LED number
    int r4 = bitRead(Led_Number, 4);   // Reads the fifth bit of the LED number

    Serial.println(Multiplex1[StringNum][0]); Serial.println(Multiplex1[StringNum][1]); Serial.println(Multiplex1[StringNum][2]); Serial.println(Multiplex1[StringNum][3]); Serial.println(Multiplex1[StringNum][4]);
    Serial.println("BLAH");

    digitalWrite(OutEnable[StringNum], LOW); // Make sure Multiplexer is turned off before switching states.

    digitalWrite(Multiplex1[StringNum][0], r0); // Switching Multiplexer to have the right state
    digitalWrite(Multiplex1[StringNum][1], r1);
    digitalWrite(Multiplex1[StringNum][2], r2);
    digitalWrite(Multiplex1[StringNum][3], r3);
    digitalWrite(Multiplex1[StringNum][4], r4);

    digitalWrite(OutEnable[StringNum], HIGH); // Turn on LED!

    delay(Time_Index); // Keep this LED on until the next thing happens
  }

}
void LED_Flash(int Led_Number, int StringNum, int Time_Index, int* OutEnable)
{
  if (Led_Number == 32)
  {
    digitalWrite(OpenLEDs[StringNum], LOW); // Turn off LED!
    delay(5);
    digitalWrite(OpenLEDs[StringNum], HIGH); // Turn on LED!
    delay(5);
    digitalWrite(OpenLEDs[StringNum], LOW); // Turn off LED!
    delay(5);
    digitalWrite(OpenLEDs[StringNum], HIGH); // Turn on LED!
    delay(5);
    digitalWrite(OpenLEDs[StringNum], LOW); // Turn on LED!
  }
  else
  {
    digitalWrite(OutEnable[StringNum], LOW); // Turn off LED!
    delay(80);
    digitalWrite(OutEnable[StringNum], HIGH); // Turn on LED!
    delay(80);
    digitalWrite(OutEnable[StringNum], LOW); // Turn off LED!
    delay(80);
    digitalWrite(OutEnable[StringNum], HIGH); // Turn on LED!
    delay(80);
    digitalWrite(OutEnable[StringNum], LOW); // Turn on LED!
    delay(Time_Index);
  }
}

