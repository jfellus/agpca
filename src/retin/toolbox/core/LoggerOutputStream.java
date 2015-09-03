package retin.toolbox.core;

import java.io.IOException;
import java.io.OutputStream;

public class LoggerOutputStream extends OutputStream {

    protected byte[] buffer;
    protected int begin,end;

    public LoggerOutputStream(int bufferSize) {
        buffer = new byte[bufferSize];
    }

    @Override
    public synchronized void write(int b) throws IOException {
        buffer[end++] = (byte)b;
        if (end >= buffer.length) {
            end = 0;
        }
        if (end == begin) {
            begin ++;
            if (begin >= buffer.length)
                begin = 0;
        }
    }

    @Override
    public synchronized String toString() {
        String res = null;
        if (end >= begin) {
            res = new String(buffer,begin,end);
        }
        else {
            res = new String(buffer,begin,buffer.length-begin)+new String(buffer,0,end);
        }
        int idx = res.indexOf('\n');
        if (idx >= 0) {
            res = res.substring(idx+1);
        }
        return res;
    }

}
