import subprocess, re, sys, os, pandas as pd, io, fcntl, time, string, datetime
from applotid import db

class AtlasI2C:
    long_timeout = 1.5 # the timeout needed to query readings and calibrations
    short_timeout = .5 # timeout for regular commands
    default_bus = 1 # the default bus for I2C on the newer Raspberry Pis, certain older boards use bus 0
    default_address = 100 # the default address for the sensor
    current_addr = default_address
    
    def __init__(self, address=default_address, bus=default_bus):
        self.file_read = io.open("/dev/i2c-"+str(bus), "rb", buffering=0)
        self.file_write = io.open("/dev/i2c-"+str(bus), "wb", buffering=0)
        self.set_i2c_address(address)

    def set_i2c_address(self, addr):
        I2C_SLAVE = 0x703
        fcntl.ioctl(self.file_read, I2C_SLAVE, addr)
        fcntl.ioctl(self.file_write, I2C_SLAVE, addr)
        self.current_addr = addr

    def write(self, cmd):
        cmd += "\00"
        self.file_write.write(bytes(cmd, 'UTF-8'))

    def read(self, num_of_bytes=32):
        res = self.file_read.read(num_of_bytes)         # read from the board
        response = list(filter(lambda x: x != '\x00', res))    # remove the null characters to get the response
        check_err = response[0]
        if check_err == 1:            
            char_list = map(lambda x: chr(x & ~0x80), list(response[1:]))
            return ''.join(char_list) 
        else:
            return "Error " + str(check_err)

    def query(self, string):
        # write a command to the board, wait the correct timeout, and read the response
        self.write(string)
        #always just do long timeout for simplicity
        time.sleep(self.long_timeout)
        return self.read()

    def close(self):
        self.file_read.close()
        self.file_write.close()
        
    def list_i2c_devices(self):
        prev_addr = self.current_addr # save the current address so we can restore it after
        i2c_devices = []
        for i in range (0,128):
            try:
                self.set_i2c_address(i)
                self.read()
                i2c_devices.append(i)
            except IOError:
                pass
        self.set_i2c_address(prev_addr) # restore the address we were using
        return i2c_devices

def dose_liquid(value):
    #Dose liquid using the dosing pump
    tentacle = AtlasI2C()
    tentacle.set_i2c_address(103)
    dispensed = tentacle.query("D,"+str(value)).rstrip("\x00")
    tentacle.close()
    return dispensed

def main():
    # Query tentacle shield: 99 = pH, 100 = E.C., 103 = Pump
    tentacle = AtlasI2C()
    # pH 
    tentacle.set_i2c_address(99)
    pH = round(float(tentacle.query("R").rstrip("\x00")),2)
    # Electrical Conductivity
    tentacle.set_i2c_address(100)
    ec_split = (tentacle.query("R").rstrip("\x00")).split(",")
    ec = ec_split[0]
    tds = ec_split[1]
    sal = ec_split[2]
    sg = ec_split[3]
    #close tentacle and save to DB
    tentacle.close()
    tentacle_dbrow = tentacleRead(datetime.datetime.now(), pH, ec, tds, sal, sg)
    db.session.add(tentacle_dbrow)
    db.session.commit()
    
if __name__ == "__main__":
    main()
   

