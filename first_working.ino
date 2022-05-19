#include <LiquidCrystal_I2C.h>
#include <SD.h>
#include <SPI.h>
#include <LedControl.h>
#include <RTClib.h>
#include <avr/pgmspace.h>

#define NEXT_BUTTON 2
#define PREV_BUTTON 3
#define MICROPHONE A0
#define MAX_SCRIPTS 8
#define MAX_SCRIPT_NAME_LEN 12
#define MATRIX_INTENSITY 2 // must be between 0 - 5
#define SAMPLES 64         // must be [16, 32, 64, ...] (64 is max for arduino UNO) 
#define SIN_ACCURACY 5     // must be between 1 - 7

struct script_t {
  File file;
  byte phaseCount;
  int phaseDelay;  
};

LiquidCrystal_I2C lcd(0x27, 16, 2);
LedControl matrix(7, 5, 6, 4);
RTC_DS3231 rtc;

bool scriptChanged = false;
byte scriptsCount = 0;
char scriptsNames[MAX_SCRIPTS][MAX_SCRIPT_NAME_LEN + 1];

byte currentScriptIndex = 0;
struct script_t currentScript;

byte currentPhase;
byte phaseBuffer[32];
int data[SAMPLES];

int lastInt = 0;

//---------------------------------  FFT stuff  ------------------------------------//
const byte isin_data[128] PROGMEM = {
  0, 1, 3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, 26, 27, 28,
  29, 31, 32, 33, 35, 36, 37, 39, 40, 41, 42, 44, 45, 46, 48, 49, 50, 52, 53, 54, 56,
  57, 59, 60, 61, 63, 64, 65, 67, 68, 70, 71, 72, 74, 75, 77, 78, 80, 81, 82, 84, 85,
  87, 88, 90, 91, 93, 94, 96, 97, 99, 100, 102, 104, 105, 107, 108, 110, 112, 113, 115,
  117, 118, 120, 122, 124, 125, 127, 129, 131, 133, 134, 136, 138, 140, 142, 144, 146,
  148, 150, 152, 155, 157, 159, 161, 164, 166, 169, 171, 174, 176, 179, 182, 185, 188,
  191, 195, 198, 202, 206, 210, 215, 221, 227, 236
};

const byte RSSdata[20] PROGMEM = {
  7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2
};

//-----------------------------FFT Function----------------------------------------------//
void Approx_FFT(int in [], int N, float Frequency) {
  int a, c1, f, o, x, data_max, data_min = 0;
  long data_avg, data_mag;
  byte scale, check = 0;

  data_max = 0;
  data_avg = 0;
  data_min = 0;

  for (int i = 0; i < 12; i++) //calculating the levels
    if ((1 << i) <= N)
      o = i;

  a = N;
  int out_r[a]; //real part of transform
  int out_im[a]; //imaginory part of transform
  
  for (int i = 0; i < a; i++) //getting min max and average for scalling
  {
    out_r[i] = 0;
    out_im[i] = 0;
    data_avg = data_avg + in [i];
    if ( in [i] > data_max) {
      data_max = in [i];
    }
    if ( in [i] < data_min) {
      data_min = in [i];
    }
  }

  data_avg = data_avg >> o;
  scale = 0;
  data_mag = data_max - data_min;

  //scalling data  from +512 to -512
  
  x = 0;
  for (int b = 0; b < o; b++) // bit reversal order stored in im_out array
  {
    c1 = (1 << b);
    f = (1 << o) / (c1 + c1);
    for (int j = 0; j < c1; j++) {
      x = x + 1;
      out_im[x] = out_im[j] + f;
    }
  }

  for (int i = 0; i < a; i++) // update input array as per bit reverse order
  {
    out_r[i] = in [out_im[i]];
    out_im[i] = 0;
  }

  int i10, i11, n1, tr, ti;
  float e;
  int c, s, temp4;
  for (int i = 0; i < o; i++) //fft
  {
    i10 = (1 << i); // overall values of sine/cosine  
    i11 = (1 << o) / (1 << (i + 1)); // loop with similar sine cosine
    e = 1024 / (1 << (i + 1)); //1024 is equivalent to 360 deg
    e = 0 - e;
    n1 = 0;

    for (int j = 0; j < i10; j++) {
      c = e * j; //c is angle as where 1024 unit is 360 deg
      while (c < 0) {
        c = c + 1024;
      }
      while (c > 1024) {
        c = c - 1024;
      }

      n1 = j;

      for (int k = 0; k < i11; k++) {
        temp4 = i10 + n1;
        if (c == 0) {
          tr = out_r[temp4];
          ti = out_im[temp4];
        } else if (c == 256) {
          tr = -out_im[temp4];
          ti = out_r[temp4];
        } else if (c == 512) {
          tr = -out_r[temp4];
          ti = -out_im[temp4];
        } else if (c == 768) {
          tr = out_im[temp4];
          ti = -out_r[temp4];
        } else if (c == 1024) {
          tr = out_r[temp4];
          ti = out_im[temp4];
        } else {
          tr = fast_cosine(out_r[temp4], c) - fast_sine(out_im[temp4], c); //the fast sine/cosine function gives direct (approx) output for A*sinx
          ti = fast_sine(out_r[temp4], c) + fast_cosine(out_im[temp4], c);
        }

        out_r[n1 + i10] = out_r[n1] - tr;
        out_r[n1] = out_r[n1] + tr;
        if (out_r[n1] > 15000 || out_r[n1] < -15000) {
          check = 1;
        } //check for int size, it can handle only +31000 to -31000,

        out_im[n1 + i10] = out_im[n1] - ti;
        out_im[n1] = out_im[n1] + ti;
        if (out_im[n1] > 15000 || out_im[n1] < -15000) {
          check = 1;
        }

        n1 = n1 + i10 + i10;
      }
    }

    if (check == 1) { // scalling the matrics if value higher than 15000 to prevent varible from overflowing
      for (int i = 0; i < a; i++) {
        out_r[i] = out_r[i] >> 1;
        out_im[i] = out_im[i] >> 1;
      }
      check = 0;
      scale = scale - 1; // tracking overall scalling of input data
    }

  }

  //---> here onward out_r contains amplitude and our_in conntains frequency (Hz) 
  int fout, fm, fstp;
  float fstep;
  fstep = Frequency / N;
  fm = 0;

  for (int i = 1; i < (1 << (o - 1)); i++) // getting amplitude from compex number
  {
    out_r[i] = fastRSS(out_r[i], out_im[i]);
    // Approx RSS function used to calculated magnitude quickly

    out_im[i] = out_im[i - 1] + fstep;
    if (fout < out_r[i]) {
      fm = i;
      fout = out_r[i];
    } 
  }

  int time = millis();
  unsigned int val;
  byte dis;
  if (out_r[1] + out_r[2] + out_r[3] < 128)
    out_r[0] = 16;
  else
    out_r[0] -= 1024;
  for (int i = 0; i < 32; i++) {
    val = abs(out_r[i]);

    if (val > 1536) {
      dis = B11111111;
    } else if (val > 1024) {
      dis = B01111111;
    } else if (val > 778) {
      dis = B00111111;
    } else if (val > 512) {
      dis = B00011111;
    } else if (val > 296) {
      dis = B00001111;
    } else if (val > 178) {
      dis = B00000111;
    } else if (val > 96) {
      dis = B00000011;
    } else if (val > 48) {
      dis = B00000001;
    } else {
      dis = B00000000;
    }

    for (byte j = 0; j < 8; j++) {
      if ((dis & (1 << j)) != 0) {
        phaseBuffer[(3 - i / 8) * 8 + 7 - j] |= (1 << (7 - i % 8));
      }
    }
  }
}

//---------------------------------fast sine/cosine---------------------------------------//

int fast_sine(int Amp, int th) {
  int temp3, m1, m2;
  byte temp1, temp2, test, quad;
  while (th > 1024) {
    th = th - 1024;
  } // here 1024 = 2*pi or 360 deg
  while (th < 0) {
    th = th + 1024;
  }
  quad = th >> 8;

  if (quad == 1) {
    th = 512 - th;
  } else if (quad == 2) {
    th = th - 512;
  } else if (quad == 3) {
    th = 1024 - th;
  }

  temp1 = 0;
  temp2 = 128; //2 multiple
  m1 = 0;
  m2 = Amp;

  temp3 = (m1 + m2) >> 1;
  Amp = temp3;
  for (int i = 0; i < SIN_ACCURACY; i++) {
    test = (temp1 + temp2) >> 1;
    temp3 = temp3 >> 1;
    if (th > pgm_read_byte(&(isin_data[test]))) {
      temp1 = test;
      Amp = Amp + temp3;
      m1 = Amp;
    } else if (th < pgm_read_byte(&(isin_data[test]))) {
      temp2 = test;
      Amp = Amp - temp3;
      m2 = Amp;
    }
  }

  if (quad == 2) {
    Amp = 0 - Amp;
  } else if (quad == 3) {
    Amp = 0 - Amp;
  }
  return (Amp);
}

int fast_cosine(int Amp, int th) {
  th = 256 - th; //cos th = sin (90-th) formula
  return (fast_sine(Amp, th));
}

//--------------------------------Fast RSS----------------------------------------//
int fastRSS(int a, int b) {
  if (a == 0 && b == 0) {
    return (0);
  }
  int min, max, temp1, temp2;
  byte clevel;
  if (a < 0) {
    a = -a;
  }
  if (b < 0) {
    b = -b;
  }
  clevel = 0;
  if (a > b) {
    max = a;
    min = b;
  } else {
    max = b;
    min = a;
  }

  if (max > (min + min + min)) {
    return max;
  } else {
    temp1 = min >> 3;
    if (temp1 == 0) {
      temp1 = 1;
    }
    temp2 = min;
    while (temp2 < max) {
      temp2 = temp2 + temp1;
      clevel = clevel + 1;
    }
    temp2 = pgm_read_byte(&(RSSdata[clevel]));
    temp1 = temp1 >> 1;
    for (int i = 0; i < temp2; i++) {
      max = max + temp1;
    }
    return (max);
  }
}
//--------------------------------------------------------------------------------//

const byte digits[60] PROGMEM = {
    0b00011000,
    0b00100100,
    0b00100100,
    0b00100100,
    0b00100100,
    0b00011000,
    
    0b00001000,
    0b00011000,
    0b00001000,
    0b00001000,
    0b00001000,
    0b00011100,
    
    0b00011000,
    0b00100100,
    0b00000100,
    0b00001000,
    0b00010000,
    0b00111100,
    
    0b00011000,
    0b00100100,
    0b00001000,
    0b00000100,
    0b00100100,
    0b00011000,
    
    0b00000100,
    0b00001100,
    0b00010100,
    0b00100100,
    0b00111110,
    0b00000100,
    
    0b00111100,
    0b00100000,
    0b00111000,
    0b00000100,
    0b00100100,
    0b00011000,
    
    0b00011000,
    0b00100000,
    0b00111000,
    0b00100100,
    0b00100100,
    0b00011000,
    
    0b00111100,
    0b00000100,
    0b00000100,
    0b00001000,
    0b00010000,
    0b00010000,
    
    0b00011000,
    0b00100100,
    0b00011000,
    0b00100100,
    0b00100100,
    0b00011000,
    
    0b00011000,
    0b00100100,
    0b00100100,
    0b00011100,
    0b00000100,
    0b00011000
};

void readScripts() {
  File scriptsDir = SD.open("scripts");
  while (true) {
    File script = scriptsDir.openNextFile();
    
    if (!script || scriptsCount == MAX_SCRIPTS)
      break;

    char *name = script.name();
    strcpy(scriptsNames[scriptsCount++], name);

    script.close();
  }
  scriptsDir.close();
}

inline int readInt(File file) {
  int acc = 0;
  int k = 1;
  byte d = file.read();

  while (d >= 48 && d <= 57) {
    acc = acc * 10 + (d - 48);
    k *= 10;
    
    d = file.read();
  }
  
  return acc;
}

inline void loadScript() {
  if (currentScript.file)
    currentScript.file.close();
  
  char path[MAX_SCRIPT_NAME_LEN + 9];
  strcpy(path, "scripts/");
  strcat(path, scriptsNames[currentScriptIndex]);

  currentScript.file = SD.open(path);
  currentScript.phaseCount = readInt(currentScript.file);
  currentScript.phaseDelay = readInt(currentScript.file);
}

inline void clearPhase() {
  for (byte i = 0; i < 32; i++)
    phaseBuffer[i] = 0;
}

inline void loadPhase() {
  char c = currentScript.file.read();
  while (c != 'I' && c != '|')
    c = currentScript.file.read();

  c = currentScript.file.read();

  if (c == 'c') {
    byte hour = rtc.now().hour();
    byte minute = rtc.now().minute();

    for (byte i = 0; i < 6; i++) {
      phaseBuffer[i + 25] = pgm_read_byte(&(digits[(hour / 10) * 6 + i])) >> 1;
      phaseBuffer[i + 17] = pgm_read_byte(&(digits[(hour % 10) * 6 + i])) << 1;
      phaseBuffer[i + 9] = pgm_read_byte(&(digits[(minute / 10) * 6 + i])) >> 1;
      phaseBuffer[i + 1] = pgm_read_byte(&(digits[(minute % 10) * 6 + i])) << 1;
    }

    for (byte i = 0; i < 4; i++) {
      phaseBuffer[i * 8] = 0;
      phaseBuffer[(i + 1) * 8 - 1] = 0;
    }
    
  } else {
    for (byte i = 0; i < 32; i++) {
      byte x = readInt(currentScript.file);
  
      for (byte j = 0; j < 8; j++) {
        if ((x & (1 << j)) != 0) {
          phaseBuffer[(3 - i / 8) * 8 + 7 - j] |= (1 << (7 - i % 8));
        }
      }
    }
  }
}

inline void showPhase() {
  for (byte i = 0; i < 32; i++)
    matrix.setRow(i / 8, i % 8, phaseBuffer[i]);
}

inline void lcdShowScript() {
  lcd.clear();
  lcd.home();
  lcd.print(currentScriptIndex + 1);
  lcd.print(". ");
  
  if (currentScriptIndex == scriptsCount)
    lcd.print("MUSIC");
  else
    lcd.print(scriptsNames[currentScriptIndex]);
}

void nextInt() {
  if (millis() - lastInt > 300) {
    currentScriptIndex += 1;
    if (currentScriptIndex > scriptsCount)
      currentScriptIndex = 0;
    scriptChanged = true;

    lastInt = millis();
  }
}

void prevInt() {
  if (millis() - lastInt > 300) {
    currentScriptIndex -= 1;
    if (currentScriptIndex > scriptsCount)
      currentScriptIndex = scriptsCount;
    scriptChanged = true;

    lastInt = millis();
  }
}


void setup() {
  SD.begin(10);
  rtc.begin();
  
  pinMode(NEXT_BUTTON, INPUT_PULLUP);
  pinMode(PREV_BUTTON, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(NEXT_BUTTON), nextInt, FALLING);
  attachInterrupt(digitalPinToInterrupt(PREV_BUTTON), prevInt, FALLING);

  lcd.init();
  lcd.clear();
  lcd.backlight();

  for (byte i = 0; i < 4; i++) {
    matrix.shutdown(i, false);
    matrix.setIntensity(i, MATRIX_INTENSITY);
    matrix.clearDisplay(i);
  }

  readScripts();
  lcdShowScript();
  
  clearPhase();
  showPhase();

  loadScript();
}

void loop() {
  if (scriptChanged) {
    scriptChanged = false;
    
    lcdShowScript();

    currentPhase = 0;
    clearPhase();
    showPhase();
    
    loadScript();
  }

  if (currentScriptIndex == scriptsCount) {
    for (byte i = 0; i < SAMPLES; i++) {
      data[i] = analogRead(MICROPHONE);
      delay(1);
    }

    clearPhase();
    Approx_FFT(data, SAMPLES, 1000.0);
    showPhase();
    
  } else if (currentScript.phaseCount > 0) {
    showPhase();

    int loadTime = millis();
    clearPhase();
    loadPhase();
    currentPhase += 1;
    if (currentPhase == currentScript.phaseCount) {
      currentPhase = 0;
      loadScript();
    }
    loadTime = millis() - loadTime;

    delay(max(currentScript.phaseDelay - loadTime, 0)); 
  } else {
    delay(300);
  }
}
