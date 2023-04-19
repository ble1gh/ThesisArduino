
#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>

// called this way, it uses the default address 0x40
Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver();
// you can also call it with a different address you want
//Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver(0x41);
// you can also call it with a different address and I2C interface
//Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver(0x40, Wire);

// Depending on your servo make, the pulse width min and max may vary, you
// want these to be as small/large as possible without hitting the hard stop
// for max range. You'll have to tweak them as necessary to match the servos you
// have!
#define SERVOMIN 100   // This is the 'minimum' pulse length count (out of 4096)
#define SERVOMAX 500   // This is the 'maximum' pulse length count (out of 4096)
#define SWEEPLEN 125   // This is the max distance from middle the servo will go
#define SERVOMID 300   // This is the 'middle' pulse length count (out of 4096)
#define USMIN 600      // This is the rounded 'minimum' microsecond length based on the minimum pulse of 150
#define USMAX 2400     // This is the rounded 'maximum' microsecond length based on the maximum pulse of 600
#define SERVO_FREQ 50  // Analog servos run at ~50 Hz updates

// our servo # counter
uint8_t servo1 = 0;
uint8_t servo2 = 1;
uint8_t servo3 = 2;
uint8_t servo4 = 3;

void setup() {
  SerialUSB.begin(115200);
  Serial.begin(115200);
  Serial.println("4 channel Servo test!");

  pwm.begin();
  /*
   * In theory the internal oscillator (clock) is 25MHz but it really isn't
   * that precise. You can 'calibrate' this by tweaking this number until
   * you get the PWM update frequency you're expecting!
   * The int.osc. for the PCA9685 chip is a range between about 23-27MHz and
   * is used for calculating things like writeMicroseconds()
   * Analog servos run at ~50 Hz updates, It is importaint to use an
   * oscilloscope in setting the int.osc frequency for the I2C PCA9685 chip.
   * 1) Attach the oscilloscope to one of the PWM signal pins and ground on
   *    the I2C PCA9685 chip you are setting the value for.
   * 2) Adjust setOscillatorFrequency() until the PWM update frequency is the
   *    expected value (50Hz for most ESCs)
   * Setting the value here is specific to each individual I2C PCA9685 chip and
   * affects the calculations for the PWM update frequency. 
   * Failure to correctly set the int.osc value will cause unexpected PWM results
   */
  pwm.setOscillatorFrequency(27000000);
  pwm.setPWMFreq(SERVO_FREQ);  // Analog servos run at ~50 Hz updates

  delay(10);
}

void loop() {

  uint16_t pos2 = 250;
  uint16_t motorvalue[4];
  uint8_t mybuffer[20] = { 0 };
  uint16_t mybufindex = 0;


  while ((SerialUSB.available() > 0) && (mybufindex < 8)) {
    uint8_t value = SerialUSB.read();
    Serial.print("I received: {");
    Serial.print(value, DEC);
    Serial.println("}");

    mybuffer[mybufindex] = value;
    mybufindex++;

    //SerialUSB.write(value);
    /*
    if ((value > SERVOMIN) && (value < SERVOMAX)) {
    motorvalue = value;
    SerialUSB.write(motorvalue);
    Serial.println("motorvalue: ");
    Serial.println(motorvalue);
    pwm.setPWM(1, 0, motorvalue);
    }    
    */
  }

  if (mybufindex == 8) {
    for (uint8_t i = 0; i < (mybufindex/2); i++) {
      motorvalue[i] = (((uint16_t)mybuffer[2*i+1]) << 8) + ((uint16_t)mybuffer[2*i]);
      Serial.print("mybuffer[2*i]: ");
      Serial.println(mybuffer[2*i]);
      Serial.print("mybuffer[2*i+1]: ");
      Serial.println(mybuffer[2*i+1]);
      Serial.print("Motor Value ");
      Serial.print(i);
      Serial.print(": ");
      Serial.println(motorvalue[i]);
    }

    for (uint8_t i = 0; i < 4; i++) {
      uint8_t highbyte = (uint8_t) ((motorvalue[i] & 0xff00)>>8);
      uint8_t lowbyte = (uint8_t) ((motorvalue[i] & 0xff));
      
      SerialUSB.write(lowbyte);
      SerialUSB.write(highbyte);
      Serial.print("lowbyte: ");
      Serial.println(lowbyte);
      Serial.print("highbyte: ");
      Serial.println(highbyte);
    }

  }

  mybufindex = 0;

  /*
  for (int i = 0; i < 4; i++) {
  	pwm.setPWM(i, 0, motorvalue);
  }
*/

  //pwm.setPWM(1, 0, pos2);
}
