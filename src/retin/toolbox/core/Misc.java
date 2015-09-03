/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package retin.toolbox.core;

import java.io.*;
import java.util.ArrayList;

/**
 *
 * @author gosselin
 */
public class Misc {

    public static String[] loadTextFile(File file) throws FileNotFoundException, IOException {
        ArrayList<String> lines = new ArrayList<String>();

        String line = null;
        BufferedReader reader = new BufferedReader(new FileReader(file));
        while ( (line=reader.readLine()) != null) {
            lines.add(line);
        }

        return lines.toArray(new String[0]);
    }


    public static void saveTextFile(File file,String content,String encoding) throws UnsupportedEncodingException, FileNotFoundException, IOException {
        Writer out = new OutputStreamWriter(new FileOutputStream(file),encoding);
        try {
            out.write(content);
        } finally {
            out.close();
        }

    }
    public static void saveTextFile(File file,String content) throws UnsupportedEncodingException, FileNotFoundException, IOException {
        Misc.saveTextFile(file,content,"UTF-8");
    }

}
